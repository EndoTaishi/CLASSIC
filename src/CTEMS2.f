      SUBROUTINE CTEMS2(FCANCMXROW,RMATCROW,ZOLNCROW,PAICROW,
     1      AILCROW,     AILCGROW,    CMASVEGCROW,  SLAICROW,
     2      AILCGSROW,   FCANCSROW,   FCANCROW,     RMATCTEMROW,
     3      CO2CONCROW,  CO2I1CGROW,  CO2I1CSROW,   CO2I2CGROW,
     4      CO2I2CSROW,  XDIFFUS,     SLAIROW,      CFLUXCGROW,
     5      CFLUXCSROW,  ANCSVEGROW,  ANCGVEGROW,   RMLCSVEGROW,
     6      RMLCGVEGROW, CANRESROW,   SDEPROW,
     7      SANDROW,     CLAYROW,     ORGMROW,
     8      ANVEGROW,    RMLVEGROW,   TCANOACCROW_M,TBARACCROW_M,
     9      UVACCROW_M,  VVACCROW_M,  MLIGHTNGGRD,  PRBFRHUCGRD,
     A      EXTNPROBGRD, STDALNGRD,   PFCANCMXROW,  NFCANCMXROW,
     B      STEMMASSROW, ROOTMASSROW, LITRMASSROW,  GLEAFMASROW,
     C      BLEAFMASROW, SOILCMASROW, AILCBROW,     FLHRLOSSROW,
     D      PANDAYSROW,  LFSTATUSROW, GRWTHEFFROW,  LYSTMMASROW,
     E      LYROTMASROW, TYMAXLAIROW, VGBIOMASROW,  GAVGLTMSROW,
     F      STMHRLOSROW, BMASVEGROW,  COLDDAYSROW,  ROTHRLOSROW,
     G      ALVSCTMROW,  ALIRCTMROW,  GAVGLAIROW,   NPPROW,
     H      NEPROW,      HETRORESROW, AUTORESROW,   SOILRESPROW,
     I      RMROW,       RGROW,       NBPROW,       LITRESROW,
     J      SOCRESROW,   GPPROW,      DSTCEMLSROW,  LITRFALLROW,
     K      HUMIFTRSROW, VEGHGHTROW,  ROOTDPTHROW,  RMLROW,
     L      RMSROW,      RMRROW,      TLTRLEAFROW,  TLTRSTEMROW,
     M      TLTRROOTROW, LEAFLITRROW, ROOTTEMPROW,  AFRLEAFROW,
     N      AFRSTEMROW,  AFRROOTROW,  WTSTATUSROW,  LTSTATUSROW,
     O      BURNFRACROW, PROBFIREROW, LUCEMCOMROW,  LUCLTRINROW,
     P      LUCSOCINROW, NPPVEGROW,   GRCLAREAROW,  DSTCEMLS3ROW,
     Q      FAREROW,     GAVGSCMSROW, TCANOACCROW_OUT,
     &      RMLVEGACCROW, RMSVEGROW,  RMRVEGROW,    RGVEGROW,
C     &      VGBIOMAS_VEGROW,GPPVEGROW,NEPVEGROW,
     &      VGBIOMAS_VEGROW,GPPVEGROW,NEPVEGROW,AILCMINROW,AILCMAXROW,
C      NEW FIRE EMISSION VARIABLES
     &      EMIT_CO2ROW,  EMIT_COROW, EMIT_CH4ROW,  EMIT_NMHCROW,
     &      EMIT_H2ROW,   EMIT_NOXROW,EMIT_N2OROW,  EMIT_PM25ROW,
     &      EMIT_TPMROW,  EMIT_TCROW, EMIT_OCROW,   EMIT_BCROW,
C  ADD FOR COMPETE
     &      LYGLFMASROW, GEREMORTROW, INTRMORTROW,  LAMBDAROW,
     &      PFTEXISTROW, BURNVEGROW,  CCROW,        MMROW,
     R      ILMOS,       JLMOS,       IWMOS,        JWMOS,
     S      NML,         NLAT,        NMOS,         ILG, 
     T      IGND,        ICAN,        ICP1,         ICC,
     U      FCANCMXGAT,  RMATCGAT,    ZOLNCGAT,     PAICGAT,
     V      AILCGAT,     AILCGGAT,    CMASVEGCGAT,  SLAICGAT,
     W      AILCGSGAT,   FCANCSGAT,   FCANCGAT,     RMATCTEMGAT,
     X      CO2CONCGAT,  CO2I1CGGAT,  CO2I1CSGAT,   CO2I2CGGAT, 
     Y      CO2I2CSGAT,  XDIFFUSGAT,  SLAIGAT,      CFLUXCGGAT, 
     Z      CFLUXCSGAT,  ANCSVEGGAT,  ANCGVEGGAT,   RMLCSVEGGAT,
     1      RMLCGVEGGAT, CANRESGAT,   SDEPGAT,
     2      SANDGAT,     CLAYGAT,     ORGMGAT,
     3      ANVEGGAT,    RMLVEGGAT,   TCANOACCGAT_M,TBARACCGAT_M,
     4      UVACCGAT_M,  VVACCGAT_M,  MLIGHTNGGAT,  PRBFRHUCGAT,
     5      EXTNPROBGAT, STDALNGAT,   PFCANCMXGAT,  NFCANCMXGAT,
     6      STEMMASSGAT, ROOTMASSGAT, LITRMASSGAT,  GLEAFMASGAT,
     7      BLEAFMASGAT, SOILCMASGAT, AILCBGAT,     FLHRLOSSGAT,
     8      PANDAYSGAT,  LFSTATUSGAT, GRWTHEFFGAT,  LYSTMMASGAT,
     9      LYROTMASGAT, TYMAXLAIGAT, VGBIOMASGAT,  GAVGLTMSGAT,
     A      STMHRLOSGAT, BMASVEGGAT,  COLDDAYSGAT,  ROTHRLOSGAT,
     B      ALVSCTMGAT,  ALIRCTMGAT,  GAVGLAIGAT,   NPPGAT,
     C      NEPGAT,      HETRORESGAT, AUTORESGAT,   SOILRESPGAT,
     D      RMGAT,       RGGAT,       NBPGAT,       LITRESGAT,
     E      SOCRESGAT,   GPPGAT,      DSTCEMLSGAT,  LITRFALLGAT,
     F      HUMIFTRSGAT, VEGHGHTGAT,  ROOTDPTHGAT,  RMLGAT,
     G      RMSGAT,      RMRGAT,      TLTRLEAFGAT,  TLTRSTEMGAT,
     H      TLTRROOTGAT, LEAFLITRGAT, ROOTTEMPGAT,  AFRLEAFGAT,
     I      AFRSTEMGAT,  AFRROOTGAT,  WTSTATUSGAT,  LTSTATUSGAT,
     J      BURNFRACGAT, PROBFIREGAT, LUCEMCOMGAT,  LUCLTRINGAT,
     K      LUCSOCINGAT, NPPVEGGAT,   GRCLAREAGAT,  DSTCEMLS3GAT,
     L      FAREGAT,     GAVGSCMSGAT, TCANOACCGAT_OUT,
     &      RMLVEGACCGAT, RMSVEGGAT,  RMRVEGGAT,    RGVEGGAT,
C     &      VGBIOMAS_VEGGAT,GPPVEGGAT,NEPVEGGAT,
     &      VGBIOMAS_VEGGAT,GPPVEGGAT,NEPVEGGAT, AILCMINGAT,AILCMAXGAT,
C      NEW FIRE EMISSION VARIABLES
     &      EMIT_CO2GAT,  EMIT_COGAT, EMIT_CH4GAT,  EMIT_NMHCGAT,
     &      EMIT_H2GAT,   EMIT_NOXGAT,EMIT_N2OGAT,  EMIT_PM25GAT,
     &      EMIT_TPMGAT,  EMIT_TCGAT, EMIT_OCGAT,   EMIT_BCGAT,
C  ADD FOR COMPETE
     &      LYGLFMASGAT, GEREMORTGAT, INTRMORTGAT,  LAMBDAGAT,
     &      PFTEXISTGAT, BURNVEGGAT,  CCGAT,        MMGAT )
C
C     * August 4/09 - Rong Li. SCATTER OPERATION ON CTEM 
C     *                         VARIABLES.
C 
      IMPLICIT NONE
C
C     * INTEGER CONSTANTS.
C
      INTEGER  NML,NLAT,NMOS,ILG,IGND,ICAN,ICP1,K,L,M
      INTEGER  ICC
C
C
C     * GATHER-SCATTER INDEX ARRAYS.
C
      INTEGER  ILMOS (ILG),  JLMOS  (ILG),  IWMOS  (ILG),  JWMOS (ILG)
C
C
      REAL  FCANCMXROW(NLAT,NMOS,ICC),RMATCROW(NLAT,NMOS,ICAN,IGND),
     1      ZOLNCROW(NLAT,NMOS,ICAN),PAICROW(NLAT,NMOS,ICAN),
     2      AILCROW(NLAT,NMOS,ICAN),AILCGROW(NLAT,NMOS,ICC),
     3      CMASVEGCROW(NLAT,NMOS,ICAN),SLAICROW(NLAT,NMOS,ICAN),
     4      AILCGSROW(NLAT,NMOS,ICC),FCANCSROW(NLAT,NMOS,ICC),
     5      FCANCROW(NLAT,NMOS,ICC),RMATCTEMROW(NLAT,NMOS,ICC,IGND),
     6      CO2CONCROW(NLAT,NMOS),CO2I1CGROW(NLAT,NMOS,ICC),
     7      CO2I1CSROW(NLAT,NMOS,ICC),CO2I2CGROW(NLAT,NMOS,ICC),
     8      CO2I2CSROW(NLAT,NMOS,ICC),XDIFFUS(NLAT),
     9      SLAIROW(NLAT,NMOS,ICC),CFLUXCGROW(NLAT,NMOS),
     A      CFLUXCSROW(NLAT,NMOS),ANCSVEGROW(NLAT,NMOS,ICC),
     B      ANCGVEGROW(NLAT,NMOS,ICC),RMLCSVEGROW(NLAT,NMOS,ICC),
     C      RMLCGVEGROW(NLAT,NMOS,ICC),CANRESROW(NLAT,NMOS),
     D      SDEPROW(NLAT,NMOS)
      REAL    SANDROW(NLAT,NMOS,IGND), CLAYROW(NLAT,NMOS,IGND), 
     1        ORGMROW(NLAT,NMOS,IGND)
      REAL  ANVEGROW(NLAT,NMOS,ICC),RMLVEGROW(NLAT,NMOS,ICC)
      REAL  TCANOACCROW_M(NLAT,NMOS),
     1      UVACCROW_M(NLAT,NMOS),VVACCROW_M(NLAT,NMOS)
      REAL  MLIGHTNGGRD(NLAT,12)  !12MONTH
      REAL  PRBFRHUCGRD(NLAT),EXTNPROBGRD(NLAT),
     1      TBARACCROW_M(NLAT,NMOS,IGND),
     2      PFCANCMXROW(NLAT,NMOS,ICC),NFCANCMXROW(NLAT,NMOS,ICC),
     3      STEMMASSROW(NLAT,NMOS,ICC),ROOTMASSROW(NLAT,NMOS,ICC),
     4      LITRMASSROW(NLAT,NMOS,ICC+1),GLEAFMASROW(NLAT,NMOS,ICC),
     5      BLEAFMASROW(NLAT,NMOS,ICC),SOILCMASROW(NLAT,NMOS,ICC+1),
     6      AILCBROW(NLAT,NMOS,ICC),FLHRLOSSROW(NLAT,NMOS,ICC)
      INTEGER   PANDAYSROW(NLAT,NMOS,ICC),LFSTATUSROW(NLAT,NMOS,ICC),
     1      STDALNGRD(NLAT),COLDDAYSROW(NLAT,NMOS,2)
      REAL  GRWTHEFFROW(NLAT,NMOS,ICC),LYSTMMASROW(NLAT,NMOS,ICC),
     9      LYROTMASROW(NLAT,NMOS,ICC),TYMAXLAIROW(NLAT,NMOS,ICC),
     A      VGBIOMASROW(NLAT,NMOS),GAVGLTMSROW(NLAT,NMOS),
     B      STMHRLOSROW(NLAT,NMOS,ICC),BMASVEGROW(NLAT,NMOS,ICC),
     C      ROTHRLOSROW(NLAT,NMOS,ICC),
     D      ALVSCTMROW(NLAT,NMOS,ICAN),ALIRCTMROW(NLAT,NMOS,ICAN),
     E      GAVGLAIROW(NLAT,NMOS)

      REAL  NPPROW(NLAT,NMOS),NEPROW(NLAT,NMOS),
     1      HETRORESROW(NLAT,NMOS),AUTORESROW(NLAT,NMOS),
     2      SOILRESPROW(NLAT,NMOS),RMROW(NLAT,NMOS),RGROW(NLAT,NMOS),
     3      NBPROW(NLAT,NMOS),LITRESROW(NLAT,NMOS),
     4      SOCRESROW(NLAT,NMOS),GPPROW(NLAT,NMOS),
     5      DSTCEMLSROW(NLAT,NMOS),LITRFALLROW(NLAT,NMOS),
     6      HUMIFTRSROW(NLAT,NMOS),VEGHGHTROW(NLAT,NMOS,ICC),
     7      ROOTDPTHROW(NLAT,NMOS,ICC),RMLROW(NLAT,NMOS),
     8      RMSROW(NLAT,NMOS),RMRROW(NLAT,NMOS),
     9      TLTRLEAFROW(NLAT,NMOS,ICC),TLTRSTEMROW(NLAT,NMOS,ICC),
     A      TLTRROOTROW(NLAT,NMOS,ICC),LEAFLITRROW(NLAT,NMOS,ICC),
     B      ROOTTEMPROW(NLAT,NMOS,ICC),AFRLEAFROW(NLAT,NMOS,ICC),
     C      AFRSTEMROW(NLAT,NMOS,ICC),AFRROOTROW(NLAT,NMOS,ICC),
     D      WTSTATUSROW(NLAT,NMOS,ICC),LTSTATUSROW(NLAT,NMOS,ICC),
     E      BURNFRACROW(NLAT,NMOS),PROBFIREROW(NLAT,NMOS),
     F      LUCEMCOMROW(NLAT,NMOS),LUCLTRINROW(NLAT,NMOS),
     G      LUCSOCINROW(NLAT,NMOS),NPPVEGROW(NLAT,NMOS,ICC),
     H      GRCLAREAROW(NLAT,NMOS),DSTCEMLS3ROW(NLAT,NMOS)
C      NEW FIRE EMISSION VARIABLES
      REAL EMIT_CO2ROW(NLAT,NMOS),  EMIT_COROW(NLAT,NMOS), 
     1     EMIT_CH4ROW(NLAT,NMOS),  EMIT_NMHCROW(NLAT,NMOS),
     2     EMIT_H2ROW(NLAT,NMOS),   EMIT_NOXROW(NLAT,NMOS),
     3     EMIT_N2OROW(NLAT,NMOS),  EMIT_PM25ROW(NLAT,NMOS),
     4     EMIT_TPMROW(NLAT,NMOS),  EMIT_TCROW(NLAT,NMOS),
     5     EMIT_OCROW(NLAT,NMOS),   EMIT_BCROW(NLAT,NMOS)

      REAL  FAREROW(NLAT,NMOS)
      REAL  GAVGSCMSROW(NLAT,NMOS)
      REAL  TCANOACCROW_OUT(NLAT,NMOS)
C
      REAL RMLVEGACCROW(NLAT,NMOS,ICC), RMSVEGROW(NLAT,NMOS,ICC),
     1      RMRVEGROW(NLAT,NMOS,ICC),    RGVEGROW(NLAT,NMOS,ICC),
     2      AILCMINROW(NLAT,NMOS,ICC),  AILCMAXROW(NLAT,NMOS,ICC)
C
      REAL VGBIOMAS_VEGROW(NLAT,NMOS,ICC)
C
      REAL GPPVEGROW(NLAT,NMOS,ICC),    NEPVEGROW(NLAT,NMOS,ICC)
C
      REAL  FCANCMXGAT(ILG,ICC),RMATCGAT(ILG,ICAN,IGND),
     1      ZOLNCGAT(ILG,ICAN),PAICGAT(ILG,ICAN),
     2      AILCGAT(ILG,ICAN),AILCGGAT(ILG,ICC),
     3      CMASVEGCGAT(ILG,ICAN),SLAICGAT(ILG,ICAN),
     4      AILCGSGAT(ILG,ICC),FCANCSGAT(ILG,ICC),
     5      FCANCGAT(ILG,ICC),RMATCTEMGAT(ILG,ICC,IGND),
     6      CO2CONCGAT(ILG),CO2I1CGGAT(ILG,ICC),
     7      CO2I1CSGAT(ILG,ICC),CO2I2CGGAT(ILG,ICC),
     8      CO2I2CSGAT(ILG,ICC),XDIFFUSGAT(ILG),
     9      SLAIGAT(ILG,ICC),CFLUXCGGAT(ILG),
     A      CFLUXCSGAT(ILG),ANCSVEGGAT(ILG,ICC),
     B      ANCGVEGGAT(ILG,ICC),RMLCSVEGGAT(ILG,ICC),
     C      RMLCGVEGGAT(ILG,ICC),CANRESGAT(ILG),
     D      SDEPGAT(ILG)
      REAL    SANDGAT(ILG,IGND), CLAYGAT(ILG,IGND), 
     1        ORGMGAT(ILG,IGND)
      REAL  ANVEGGAT(ILG,ICC),RMLVEGGAT(ILG,ICC)
      REAL  TCANOACCGAT_M(ILG),
     1      UVACCGAT_M(ILG),VVACCGAT_M(ILG)
      REAL  MLIGHTNGGAT(ILG,12)  !12MONTH
      REAL  PRBFRHUCGAT(ILG),EXTNPROBGAT(ILG),
     1      TBARACCGAT_M(ILG,IGND),
     2      PFCANCMXGAT(ILG,ICC),NFCANCMXGAT(ILG,ICC),
     3      STEMMASSGAT(ILG,ICC),  ROOTMASSGAT(ILG,ICC),  
     4      LITRMASSGAT(ILG,ICC+1),GLEAFMASGAT(ILG,ICC),
     5      BLEAFMASGAT(ILG,ICC),SOILCMASGAT(ILG,ICC+1),
     6      AILCBGAT(ILG,ICC),FLHRLOSSGAT(ILG,ICC)
      INTEGER     PANDAYSGAT(ILG,ICC),LFSTATUSGAT(ILG,ICC),
     1      STDALNGAT(ILG),COLDDAYSGAT(ILG,2)
      REAL  GRWTHEFFGAT(ILG,ICC),LYSTMMASGAT(ILG,ICC),
     9      LYROTMASGAT(ILG,ICC),TYMAXLAIGAT(ILG,ICC),
     A      VGBIOMASGAT(ILG),GAVGLTMSGAT(ILG),
     B      STMHRLOSGAT(ILG,ICC),BMASVEGGAT(ILG,ICC),
     C      ROTHRLOSGAT(ILG,ICC),
     D      ALVSCTMGAT(ILG,ICAN),ALIRCTMGAT(ILG,ICAN),
     E      GAVGLAIGAT(ILG)
C
      REAL  NPPGAT(ILG),NEPGAT(ILG),
     1      HETRORESGAT(ILG),AUTORESGAT(ILG),
     2      SOILRESPGAT(ILG),RMGAT(ILG),RGGAT(ILG),
     3      NBPGAT(ILG),LITRESGAT(ILG),
     4      SOCRESGAT(ILG),GPPGAT(ILG),
     5      DSTCEMLSGAT(ILG),LITRFALLGAT(ILG),
     6      HUMIFTRSGAT(ILG),VEGHGHTGAT(ILG,ICC),
     7      ROOTDPTHGAT(ILG,ICC),RMLGAT(ILG),
     8      RMSGAT(ILG),RMRGAT(ILG),
     9      TLTRLEAFGAT(ILG,ICC),TLTRSTEMGAT(ILG,ICC),
     A      TLTRROOTGAT(ILG,ICC),LEAFLITRGAT(ILG,ICC),
     B      ROOTTEMPGAT(ILG,ICC),AFRLEAFGAT(ILG,ICC),
     C      AFRSTEMGAT(ILG,ICC),AFRROOTGAT(ILG,ICC),
     D      WTSTATUSGAT(ILG,ICC),LTSTATUSGAT(ILG,ICC),
     E      BURNFRACGAT(ILG),PROBFIREGAT(ILG),
     F      LUCEMCOMGAT(ILG),LUCLTRINGAT(ILG),
     G      LUCSOCINGAT(ILG),NPPVEGGAT(ILG,ICC),
     H      GRCLAREAGAT(ILG),DSTCEMLS3GAT(ILG)
C      NEW FIRE EMISSION VARIABLES
       REAL EMIT_CO2GAT(ILG),  EMIT_COGAT(ILG), 
     1      EMIT_CH4GAT(ILG),  EMIT_NMHCGAT(ILG),
     2      EMIT_H2GAT(ILG),   EMIT_NOXGAT(ILG),
     3      EMIT_N2OGAT(ILG),  EMIT_PM25GAT(ILG),
     4      EMIT_TPMGAT(ILG),  EMIT_TCGAT(ILG),
     5      EMIT_OCGAT(ILG),   EMIT_BCGAT(ILG)

      REAL  FAREGAT(ILG) 
      REAL  GAVGSCMSGAT(ILG)  
      REAL  TCANOACCGAT_OUT(ILG)
C
      REAL RMLVEGACCGAT(ILG,ICC),     RMSVEGGAT(ILG,ICC),
     1      RMRVEGGAT(ILG,ICC),       RGVEGGAT(ILG,ICC),
     2      AILCMINGAT(ILG,ICC),     AILCMAXGAT(ILG,ICC)
C
      REAL VGBIOMAS_VEGGAT(ILG,ICC)
C
      REAL GPPVEGGAT(ILG,ICC),        NEPVEGGAT(ILG,ICC)

C  ADD FOR COMPETE
       REAL LYGLFMASROW(NLAT,NMOS,ICC), LYGLFMASGAT(ILG,ICC),
     1      GEREMORTROW(NLAT,NMOS,ICC), GEREMORTGAT(ILG,ICC),
     2      INTRMORTROW(NLAT,NMOS,ICC), INTRMORTGAT(ILG,ICC),
     3      LAMBDAROW(NLAT,NMOS,ICC), LAMBDAGAT(ILG,ICC),
     4      PFTEXISTROW(NLAT,NMOS,ICC), PFTEXISTGAT(ILG,ICC),
     5      BURNVEGROW(NLAT,NMOS,ICC), BURNVEGGAT(ILG,ICC),
     6      CCROW(NLAT,NMOS,ICC), CCGAT(ILG,ICC),
     7      MMROW(NLAT,NMOS,ICC), MMGAT(ILG,ICC)
C
C----------------------------------------------------------------------
      DO 100 K=1,NML
	  SDEPROW(ILMOS(K),JLMOS(K))=SDEPGAT(K)
	  CO2CONCROW(ILMOS(K),JLMOS(K))=CO2CONCGAT(K)
	  CFLUXCGROW(ILMOS(K),JLMOS(K))=CFLUXCGGAT(K)
          CFLUXCSROW(ILMOS(K),JLMOS(K))=CFLUXCSGAT(K)
          CANRESROW(ILMOS(K),JLMOS(K))=CANRESGAT(K)
          XDIFFUS(ILMOS(K))=XDIFFUSGAT(K) 
          PRBFRHUCGRD(ILMOS(K))=PRBFRHUCGAT(K) 
          EXTNPROBGRD(ILMOS(K))=EXTNPROBGAT(K) 
          STDALNGRD(ILMOS(K))=STDALNGAT(K) 
          TCANOACCROW_M(ILMOS(K),JLMOS(K))=TCANOACCGAT_M(K) 
          UVACCROW_M(ILMOS(K),JLMOS(K))=UVACCGAT_M(K) 
          VVACCROW_M(ILMOS(K),JLMOS(K))=VVACCGAT_M(K) 
          VGBIOMASROW(ILMOS(K),JLMOS(K))=VGBIOMASGAT(K)
          GAVGLTMSROW(ILMOS(K),JLMOS(K))=GAVGLTMSGAT(K)
          GAVGLAIROW(ILMOS(K),JLMOS(K))=GAVGLAIGAT(K)
          NPPROW(ILMOS(K),JLMOS(K))=NPPGAT(K)
          NEPROW(ILMOS(K),JLMOS(K))=NEPGAT(K)
          HETRORESROW(ILMOS(K),JLMOS(K))=HETRORESGAT(K)
          AUTORESROW(ILMOS(K),JLMOS(K))=AUTORESGAT(K)
          SOILRESPROW(ILMOS(K),JLMOS(K))=SOILRESPGAT(K)
          RMROW(ILMOS(K),JLMOS(K))=RMGAT(K)
          RGROW(ILMOS(K),JLMOS(K))=RGGAT(K)
          NBPROW(ILMOS(K),JLMOS(K))=NBPGAT(K)
          LITRESROW(ILMOS(K),JLMOS(K))=LITRESGAT(K)
          SOCRESROW(ILMOS(K),JLMOS(K))=SOCRESGAT(K)
          GPPROW(ILMOS(K),JLMOS(K))=GPPGAT(K)
          DSTCEMLSROW(ILMOS(K),JLMOS(K))=DSTCEMLSGAT(K)
          LITRFALLROW(ILMOS(K),JLMOS(K))=LITRFALLGAT(K)
          HUMIFTRSROW(ILMOS(K),JLMOS(K))=HUMIFTRSGAT(K)
          RMLROW(ILMOS(K),JLMOS(K))=RMLGAT(K)
          RMSROW(ILMOS(K),JLMOS(K))=RMSGAT(K)
          RMRROW(ILMOS(K),JLMOS(K))=RMRGAT(K)
          BURNFRACROW(ILMOS(K),JLMOS(K))=BURNFRACGAT(K)
          PROBFIREROW(ILMOS(K),JLMOS(K))=PROBFIREGAT(K)
          LUCEMCOMROW(ILMOS(K),JLMOS(K))=LUCEMCOMGAT(K)
          LUCLTRINROW(ILMOS(K),JLMOS(K))=LUCLTRINGAT(K)
          LUCSOCINROW(ILMOS(K),JLMOS(K))=LUCSOCINGAT(K)
          GRCLAREAROW(ILMOS(K),JLMOS(K))=GRCLAREAGAT(K)
          DSTCEMLS3ROW(ILMOS(K),JLMOS(K))=DSTCEMLS3GAT(K)
          FAREROW(ILMOS(K),JLMOS(K))=FAREGAT(K)
          GAVGSCMSROW(ILMOS(K),JLMOS(K))=GAVGSCMSGAT(k)
          TCANOACCROW_OUT(ILMOS(K),JLMOS(K))=TCANOACCGAT_OUT(K)
C      NEW FIRE EMISSION VARIABLES
          EMIT_CO2ROW(ILMOS(K),JLMOS(K))=EMIT_CO2GAT(K)
          EMIT_COROW(ILMOS(K),JLMOS(K))=EMIT_COGAT(K) 
          EMIT_CH4ROW(ILMOS(K),JLMOS(K))=EMIT_CH4GAT(K)
          EMIT_NMHCROW(ILMOS(K),JLMOS(K))=EMIT_NMHCGAT(K)
          EMIT_H2ROW(ILMOS(K),JLMOS(K))=EMIT_H2GAT(K)
          EMIT_NOXROW(ILMOS(K),JLMOS(K))=EMIT_NOXGAT(K)
          EMIT_N2OROW(ILMOS(K),JLMOS(K))=EMIT_N2OGAT(K)
          EMIT_PM25ROW(ILMOS(K),JLMOS(K))=EMIT_PM25GAT(K)
          EMIT_TPMROW(ILMOS(K),JLMOS(K))=EMIT_TPMGAT(K)
          EMIT_TCROW(ILMOS(K),JLMOS(K))=EMIT_TCGAT(K)
          EMIT_OCROW(ILMOS(K),JLMOS(K))=EMIT_OCGAT(K)
          EMIT_BCROW(ILMOS(K),JLMOS(K))= EMIT_BCGAT(K)
100   CONTINUE
C
      DO 101 L=1,ICC
      DO 101 K=1,NML
	  AILCGROW(ILMOS(K),JLMOS(K),L)=AILCGGAT(K,L)   
          AILCGSROW(ILMOS(K),JLMOS(K),L)=AILCGSGAT(K,L)
          CO2I1CGROW(ILMOS(K),JLMOS(K),L)=CO2I1CGGAT(K,L) 
          CO2I1CSROW(ILMOS(K),JLMOS(K),L)=CO2I1CSGAT(K,L) 
          CO2I2CGROW(ILMOS(K),JLMOS(K),L)=CO2I2CGGAT(K,L) 
          CO2I2CSROW(ILMOS(K),JLMOS(K),L)=CO2I2CSGAT(K,L) 
          SLAIROW(ILMOS(K),JLMOS(K),L)=SLAIGAT(K,L) 
          ANVEGROW(ILMOS(K),JLMOS(K),L)=ANVEGGAT(K,L) 
          RMLVEGROW(ILMOS(K),JLMOS(K),L)=RMLVEGGAT(K,L) 

C         Vivek ---\\
          PFCANCMXROW(ILMOS(K),JLMOS(K),L) = PFCANCMXGAT(K,L)
          NFCANCMXROW(ILMOS(K),JLMOS(K),L) = NFCANCMXGAT(K,L)
          FCANCMXROW(ILMOS(K),JLMOS(K),L) = FCANCMXGAT(K,L)
C         Vivek ---//

	      STEMMASSROW(ILMOS(K),JLMOS(K),L)=STEMMASSGAT(K,L)
          ROOTMASSROW(ILMOS(K),JLMOS(K),L)=ROOTMASSGAT(K,L)
          GLEAFMASROW(ILMOS(K),JLMOS(K),L)=GLEAFMASGAT(K,L)
      	  BLEAFMASROW(ILMOS(K),JLMOS(K),L)=BLEAFMASGAT(K,L)
          AILCBROW(ILMOS(K),JLMOS(K),L)=AILCBGAT(K,L)   
          FLHRLOSSROW(ILMOS(K),JLMOS(K),L)=FLHRLOSSGAT(K,L)
          PANDAYSROW(ILMOS(K),JLMOS(K),L)=PANDAYSGAT(K,L)
          LFSTATUSROW(ILMOS(K),JLMOS(K),L)=LFSTATUSGAT(K,L)
          GRWTHEFFROW(ILMOS(K),JLMOS(K),L)=GRWTHEFFGAT(K,L)
          LYSTMMASROW(ILMOS(K),JLMOS(K),L)=LYSTMMASGAT(K,L)
          LYROTMASROW(ILMOS(K),JLMOS(K),L)=LYROTMASGAT(K,L)
          TYMAXLAIROW(ILMOS(K),JLMOS(K),L)=TYMAXLAIGAT(K,L)
          STMHRLOSROW(ILMOS(K),JLMOS(K),L)=STMHRLOSGAT(K,L)
          BMASVEGROW(ILMOS(K),JLMOS(K),L)=BMASVEGGAT(K,L)
          ROTHRLOSROW(ILMOS(K),JLMOS(K),L)=ROTHRLOSGAT(K,L)
          VEGHGHTROW(ILMOS(K),JLMOS(K),L)=VEGHGHTGAT(K,L)
          ROOTDPTHROW(ILMOS(K),JLMOS(K),L)=ROOTDPTHGAT(K,L)
          TLTRLEAFROW(ILMOS(K),JLMOS(K),L)=TLTRLEAFGAT(K,L)
          TLTRSTEMROW(ILMOS(K),JLMOS(K),L)=TLTRSTEMGAT(K,L)
          TLTRROOTROW(ILMOS(K),JLMOS(K),L)=TLTRROOTGAT(K,L)
          LEAFLITRROW(ILMOS(K),JLMOS(K),L)=LEAFLITRGAT(K,L)
          ROOTTEMPROW(ILMOS(K),JLMOS(K),L)=ROOTTEMPGAT(K,L)
          AFRLEAFROW(ILMOS(K),JLMOS(K),L)=AFRLEAFGAT(K,L)
          AFRSTEMROW(ILMOS(K),JLMOS(K),L)=AFRSTEMGAT(K,L)
          AFRROOTROW(ILMOS(K),JLMOS(K),L)=AFRROOTGAT(K,L)
          WTSTATUSROW(ILMOS(K),JLMOS(K),L)=WTSTATUSGAT(K,L)
          LTSTATUSROW(ILMOS(K),JLMOS(K),L)=LTSTATUSGAT(K,L)
          NPPVEGROW(ILMOS(K),JLMOS(K),L)=NPPVEGGAT(K,L)
          RMLVEGACCROW(ILMOS(K),JLMOS(K),L)=RMLVEGACCGAT(K,L)
          RMSVEGROW(ILMOS(K),JLMOS(K),L)=RMSVEGGAT(K,L)
          RMRVEGROW(ILMOS(K),JLMOS(K),L)=RMRVEGGAT(K,L)
          RGVEGROW(ILMOS(K),JLMOS(K),L)=RGVEGGAT(K,L)
          VGBIOMAS_VEGROW(ILMOS(K),JLMOS(K),L)=VGBIOMAS_VEGGAT(K,L)
          GPPVEGROW(ILMOS(K),JLMOS(K),L)=GPPVEGGAT(K,L)
          NEPVEGROW(ILMOS(K),JLMOS(K),L)=NEPVEGGAT(K,L)
C FLAG these were missing previously JM  
          AILCMINROW(ILMOS(K),JLMOS(K),L)=AILCMINGAT(K,L)
          AILCMAXROW(ILMOS(K),JLMOS(K),L)=AILCMAXGAT(K,L)
C
C  ADD FOR COMPETE
          LYGLFMASROW(ILMOS(K),JLMOS(K),L)=LYGLFMASGAT(K,L)
          GEREMORTROW(ILMOS(K),JLMOS(K),L)=GEREMORTGAT(K,L)
          INTRMORTROW(ILMOS(K),JLMOS(K),L)=INTRMORTGAT(K,L)
          LAMBDAROW(ILMOS(K),JLMOS(K),L)=LAMBDAGAT(K,L)
          PFTEXISTROW(ILMOS(K),JLMOS(K),L)=PFTEXISTGAT(K,L)
          BURNVEGROW(ILMOS(K),JLMOS(K),L)=BURNVEGGAT(K,L)
          CCROW(ILMOS(K),JLMOS(K),L)=CCGAT(K,L)
          MMROW(ILMOS(K),JLMOS(K),L)=MMGAT(K,L)
C
101   CONTINUE
C
      DO 102 L=1,ICC+1
      DO 102 K=1,NML
	  LITRMASSROW(ILMOS(K),JLMOS(K),L)=LITRMASSGAT(k,L)
	  SOILCMASROW(ILMOS(K),JLMOS(K),L)=SOILCMASGAT(k,L)
102   CONTINUE
C
      DO 105 L=1,12     !12 MONTHS
      DO 105 K=1,NML
          MLIGHTNGGRD(ILMOS(K),L)=MLIGHTNGGAT(K,L)
105   CONTINUE
C
      DO 106 L=1,2     !2 PFTs (NDL DCD & CROPs)
      DO 106 K=1,NML
          COLDDAYSROW(ILMOS(K),JLMOS(K),L)=COLDDAYSGAT(K,L)
106   CONTINUE
C
      DO 201 L=1,ICAN
      DO 201 K=1,NML
          AILCROW(ILMOS(K),JLMOS(K),L)=AILCGAT(K,L)
          ZOLNCROW(ILMOS(K),JLMOS(K),L)=ZOLNCGAT(K,L)
	  CMASVEGCROW(ILMOS(K),JLMOS(K),L)=CMASVEGCGAT(K,L)
          PAICROW(ILMOS(K),JLMOS(K),L)=PAICGAT(K,L)
          SLAICROW(ILMOS(K),JLMOS(K),L)=SLAICGAT(K,L)
          ALVSCTMROW(ILMOS(K),JLMOS(K),L)=ALVSCTMGAT(K,L)
          ALIRCTMROW(ILMOS(K),JLMOS(K),L)=ALIRCTMGAT(K,L)
201   CONTINUE
C
      DO 250 L=1,IGND
      DO 250 K=1,NML
          SANDROW(ILMOS(K),JLMOS(K),L)=SANDGAT(K,L)
          CLAYROW(ILMOS(K),JLMOS(K),L)=CLAYGAT(K,L)
          ORGMROW(ILMOS(K),JLMOS(K),L)=ORGMGAT(K,L)
          TBARACCROW_M(ILMOS(K),JLMOS(K),L)=TBARACCGAT_M(K,L)
250   CONTINUE
C
      DO 280 L=1,ICC
      DO 280 M=1,IGND
      DO 280 K=1,NML
	  RMATCTEMROW(ILMOS(K),JLMOS(K),L,M)=RMATCTEMGAT(K,L,M)
280   CONTINUE
C
      DO 290 L=1,ICAN
      DO 290 M=1,IGND
      DO 290 K=1,NML
	  RMATCROW(ILMOS(K),JLMOS(K),L,M)=RMATCGAT(K,L,M)
290   CONTINUE
C
      RETURN
      END
