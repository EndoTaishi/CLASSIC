      SUBROUTINE CTEMG1(  GLEAFMASGAT,  BLEAFMASGAT,  STEMMASSGAT,
     A                    ROOTMASSGAT,   FCANCMXGAT,      ZBTWGAT,
     B                        DLZWGAT,      SDEPGAT,     AILCGGAT,
     C           AILCBGAT,    AILCGAT,     ZOLNCGAT,     RMATCGAT,
     D        RMATCTEMGAT,    SLAIGAT,   BMASVEGGAT,  CMASVEGCGAT,
     E         VEGHGHTGAT,ROOTDPTHGAT,   ALVSCTMGAT,   ALIRCTMGAT,
     F            PAICGAT,   SLAICGAT,     SANDGAT,       CLAYGAT,
     G              ILMOS,      JLMOS,       IWMOS,         JWMOS,
     H         NML,  NLAT,  NMOS, ILG,  IGND, ICAN,  ICP1,    ICC,
     I        GLEAFMASROW,BLEAFMASROW, STEMMASSROW,   ROOTMASSROW,
     J         FCANCMXROW,    ZBTWROW,     DLZWROW,       SDEPROW,
     K           AILCGROW,   AILCBROW,     AILCROW,      ZOLNCROW,
     L           RMATCROW,RMATCTEMROW,     SLAIROW,    BMASVEGROW,
     M        CMASVEGCROW, VEGHGHTROW, ROOTDPTHROW,    ALVSCTMROW,
     N         ALIRCTMROW,    PAICROW,    SLAICROW,       SANDROW,
     O            CLAYROW)
C
C     July 5 2009   GATHER OPERATION ON CTEM VARIABLES FOR CONSISTENCY
C                   WITH CLASS' TILED VERSION
C     RONG LI
C
      IMPLICIT NONE
C
C     INTEGER CONSTANTS.
C
      INTEGER  NML, NLAT, NMOS, ILG, IGND, ICAN, ICP1, K, L, M
      INTEGER  ICC
C
C
C     GATHER-SCATTER INDEX ARRAYS.
C
      INTEGER  ILMOS (ILG),  JLMOS  (ILG),  IWMOS  (ILG),  JWMOS (ILG)
C
C
      REAL  GLEAFMASGAT(ILG,ICC),  BLEAFMASGAT(ILG,ICC),  
     1      STEMMASSGAT(ILG,ICC),  ROOTMASSGAT(ILG,ICC),  
     2       FCANCMXGAT(ILG,ICC),     ZBTWGAT(ILG,IGND), 
     3         DLZWGAT(ILG,IGND),          SDEPGAT(ILG),      
     4         AILCGGAT(ILG,ICC),     AILCBGAT(ILG,ICC),
     5         AILCGAT(ILG,ICAN),     ZOLNCGAT(ILG,ICAN),
     6   RMATCGAT(ILG,ICAN,IGND),
     7 RMATCTEMGAT(ILG,ICC,IGND),
     8          SLAIGAT(ILG,ICC),    BMASVEGGAT(ILG,ICC),
     9     CMASVEGCGAT(ILG,ICAN),    VEGHGHTGAT(ILG,ICC),
     A      ROOTDPTHGAT(ILG,ICC),
     B      ALVSCTMGAT(ILG,ICAN),   ALIRCTMGAT(ILG,ICAN),
     C         PAICGAT(ILG,ICAN),     SLAICGAT(ILG,ICAN),
     D         SANDGAT(ILG,IGND),      CLAYGAT(ILG,IGND)
C
      REAL  GLEAFMASROW(NLAT,NMOS,ICC),     BLEAFMASROW(NLAT,NMOS,ICC),
     1      STEMMASSROW(NLAT,NMOS,ICC),     ROOTMASSROW(NLAT,NMOS,ICC),
     2       FCANCMXROW(NLAT,NMOS,ICC),        ZBTWROW(NLAT,NMOS,IGND),
     3         DLZWROW(NLAT,NMOS,IGND),             SDEPROW(NLAT,NMOS),
     4         AILCGROW(NLAT,NMOS,ICC),        AILCBROW(NLAT,NMOS,ICC),
     5         AILCROW(NLAT,NMOS,ICAN),       ZOLNCROW(NLAT,NMOS,ICAN),
     6   RMATCROW(NLAT,NMOS,ICAN,IGND),
     7 RMATCTEMROW(NLAT,NMOS,ICC,IGND),
     8          SLAIROW(NLAT,NMOS,ICC),      BMASVEGROW(NLAT,NMOS,ICC),
     9     CMASVEGCROW(NLAT,NMOS,ICAN),      VEGHGHTROW(NLAT,NMOS,ICC),
     A      ROOTDPTHROW(NLAT,NMOS,ICC),
     B      ALVSCTMROW(NLAT,NMOS,ICAN),     ALIRCTMROW(NLAT,NMOS,ICAN),
     C         PAICROW(NLAT,NMOS,ICAN),       SLAICROW(NLAT,NMOS,ICAN),
     D         SANDROW(NLAT,NMOS,IGND),        CLAYROW(NLAT,NMOS,IGND)
C
C----------------------------------------------------------------------
      DO 100 K=1,NML
          SDEPGAT(K)=SDEPROW(ILMOS(K),JLMOS(K))
100   CONTINUE
C
      DO 101 L=1,ICC
       DO 101 K=1,NML
          GLEAFMASGAT(K,L) = GLEAFMASROW(ILMOS(K),JLMOS(K),L)
          BLEAFMASGAT(K,L) = BLEAFMASROW(ILMOS(K),JLMOS(K),L)
          STEMMASSGAT(K,L) = STEMMASSROW(ILMOS(K),JLMOS(K),L)
          ROOTMASSGAT(K,L) = ROOTMASSROW(ILMOS(K),JLMOS(K),L)
          FCANCMXGAT(K,L)  = FCANCMXROW(ILMOS(K),JLMOS(K),L)
          AILCGGAT(K,L)    = AILCGROW(ILMOS(K),JLMOS(K),L)
          AILCBGAT(K,L)    = AILCBROW(ILMOS(K),JLMOS(K),L)
          SLAIGAT(K,L)     = SLAIROW(ILMOS(K),JLMOS(K),L)
          BMASVEGGAT(K,L)  = BMASVEGROW(ILMOS(K),JLMOS(K),L)
          VEGHGHTGAT(K,L)  = VEGHGHTROW(ILMOS(K),JLMOS(K),L)
          ROOTDPTHGAT(K,L) = ROOTDPTHROW(ILMOS(K),JLMOS(K),L)
101   CONTINUE
C
      DO 201 L=1,ICAN
       DO 201 K=1,NML
          AILCGAT(K,L)     = AILCROW(ILMOS(K),JLMOS(K),L)
          ZOLNCGAT(K,L)    = ZOLNCROW(ILMOS(K),JLMOS(K),L)
          CMASVEGCGAT(K,L) = CMASVEGCROW(ILMOS(K),JLMOS(K),L)
          ALVSCTMGAT(K,L)  = ALVSCTMROW(ILMOS(K),JLMOS(K),L)
          ALIRCTMGAT(K,L)  = ALIRCTMROW(ILMOS(K),JLMOS(K),L)
          PAICGAT(K,L)     = PAICROW(ILMOS(K),JLMOS(K),L)
          SLAICGAT(K,L)    = SLAICROW(ILMOS(K),JLMOS(K),L)
201   CONTINUE
C
      DO 250 L=1,IGND
       DO 250 K=1,NML
          ZBTWGAT(K,L) = ZBTWROW(ILMOS(K),JLMOS(K),L)
          DLZWGAT(K,L) = DLZWROW(ILMOS(K),JLMOS(K),L)
          SANDGAT(K,L) = SANDROW(ILMOS(K),JLMOS(K),L)
          CLAYGAT(K,L) = CLAYROW(ILMOS(K),JLMOS(K),L)
250   CONTINUE
C
      DO 280 L=1,ICC
       DO 280 M=1,IGND
        DO 280 K=1,NML
          RMATCTEMGAT(K,L,M) = RMATCTEMROW(ILMOS(K),JLMOS(K),L,M)
280   CONTINUE
C
      DO 290 L=1,ICAN
       DO 290 M=1,IGND
        DO 290 K=1,NML
           RMATCGAT(K,L,M)=RMATCROW(ILMOS(K),JLMOS(K),L,M)
290   CONTINUE
C
      RETURN
      END
