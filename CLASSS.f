      SUBROUTINE CLASSS (TBARROT,THLQROT,THICROT,TSFSROT,TPNDROT,       
     1                   ZPNDROT,TBASROT,ALBSROT,TSNOROT,RHOSROT,       
     2                   SNOROT, GTROT,  TCANROT,RCANROT,SCANROT,       
     3                  GROROT, CMAIROT,TACROT, QACROT, WSNOROT,
     4                   REFROT, BCSNROT,EMISROT,SALBROT,CSALROT,       
     5                   ILMOS,JLMOS,                                   
     6                   NML,NL,NT,NM,ILG,IG,IC,ICP1,NBS,               
     7                   TBARGAT,THLQGAT,THICGAT,TSFSGAT,TPNDGAT,       
     8                   ZPNDGAT,TBASGAT,ALBSGAT,TSNOGAT,RHOSGAT,       
     9                   SNOGAT, GTGAT,  TCANGAT,RCANGAT,SCANGAT,       
     A                   GROGAT, CMAIGAT,TACGAT, QACGAT, WSNOGAT,       
     B                   REFGAT, BCSNGAT,EMISGAT,SALBGAT,CSALGAT)       
C                                                                       
C     * Jun 20, 2014 - M.Lazare. New version for gcm18, called          
C     *                          by new "sfcproc2":                     
C     *                          - Adds SALBGAT/SALBROT and             
C     *                            CSALGAT,CSALROT (need to pass        
C     *                            NBS as well).                        
C     *                          - Adds EMISGAT/EMISROT.                
C     *                          - Adds GTGAT/GTROT.                    
C     *                          - Adds NT (NTLD in sfcproc2            
C     *                            call) to dimension land-only         
C     *                            ROT fields, consistent with          
C     *                            new comrow12.                        
C     *                          - Unused IWMOS,JWMOS removed.          
C     * Jun 12, 2013 - M.Lazare. Previous version for gcm17,            
C     *                          called by "sfcproc".                   
C     *                          CLASS scatter routine called by        
C     *                          "sfcproc" in new version gcm17.        
C     * NOTE: This contains the following changes compared to the       
C     *       working temporary version used in conjunction with        
C     *       updates to gcm16 (ie not official):                       
C     *         1) {REF,BCSN} added for Maryam's new code.              
C     *         2) GFLX removed.                                        
C
C     * OCT 25/11 - M.LAZARE.   REMOVE OPERATIONS ON INTERNAL
C     *                         ROT ARRAYS (NOW DONE DIRECTLY
C     *                         GAT->ROW IN SFCPROC).
C     * OCT 07/11 - M.LAZARE.   REMOVE TSF.
C     * OCT 05/11 - M.LAZARE.   ADD SFCH.
C     * OCT 04/11 - M.LAZARE.   REMOVE ITCT.
C     * MAR 23/06 - D.VERSEGHY. ADD WSNO,FSNO.
C     * MAR 18/05 - D.VERSEGHY. ADDITIONAL VARIABLES.
C     * FEB 18/05 - D.VERSEGHY. ADD "TSFS" VARIABLES.
C     * AUG 05/04 - D.VERSEGHY. ADD NEW DIAGNOSTIC VARIABLES
C     *                         ILMO, UE AND HBL.
C     * AUG 15/02 - D.VERSEGHY. SCATTER OPERATION ON CLASS 
C     *                         VARIABLES.
C
      IMPLICIT NONE
C
C     * INTEGER CONSTANTS.
C
      INTEGER NML,NL,NT,NM,ILG,IG,IC,ICP1,NBS,K,L,M                     
C
C     * LAND SURFACE PROGNOSTIC VARIABLES.
C
      REAL    SALBROT(NL,NM,NBS),CSALROT(NL,NM,NBS)                     
      REAL    TBARROT(NL,NT,IG), THLQROT(NL,NT,IG), THICROT(NL,NT,IG)   
      REAL    TSFSROT(NL,NT,4)                                          
      REAL    TPNDROT(NL,NT),    ZPNDROT(NL,NT),    TBASROT(NL,NT),     
     1        ALBSROT(NL,NM),    TSNOROT(NL,NM),    RHOSROT(NL,NM),   
     2        SNOROT (NL,NM),    GTROT  (NL,NM),    TCANROT(NL,NT),     
     3        RCANROT(NL,NT),    SCANROT(NL,NT),    GROROT (NL,NT),     
     4        TACROT (NL,NT),    QACROT (NL,NT),    WSNOROT(NL,NT),     
     5        CMAIROT(NL,NT),    REFROT (NL,NM),    BCSNROT(NL,NM),     
     6        EMISROT(NL,NM)                                            
C
      REAL    SALBGAT(ILG,NBS),  CSALGAT(ILG,NBS)                       
      REAL    TBARGAT(ILG,IG),   THLQGAT(ILG,IG),   THICGAT(ILG,IG)     

      REAL    TSFSGAT(ILG,4)

      REAL    TPNDGAT(ILG),      ZPNDGAT(ILG),      TBASGAT(ILG),   
     1        ALBSGAT(ILG),      TSNOGAT(ILG),      RHOSGAT(ILG),   
     2        SNOGAT (ILG),      GTGAT  (ILG),      TCANGAT(ILG),       
     3        RCANGAT(ILG),      SCANGAT(ILG),      GROGAT (ILG),       
     4        TACGAT (ILG),      QACGAT (ILG),      WSNOGAT(ILG),       
     5        CMAIGAT(ILG),      REFGAT (ILG),      BCSNGAT(ILG),       
     6        EMISGAT(ILG)                                              
C
C     * GATHER-SCATTER INDEX ARRAYS.
C
      INTEGER  ILMOS (ILG),  JLMOS  (ILG)                               
C----------------------------------------------------------------------
      DO 100 K=1,NML
          TPNDROT(ILMOS(K),JLMOS(K))=TPNDGAT(K)  
          ZPNDROT(ILMOS(K),JLMOS(K))=ZPNDGAT(K)  
          TBASROT(ILMOS(K),JLMOS(K))=TBASGAT(K)  
          ALBSROT(ILMOS(K),JLMOS(K))=ALBSGAT(K)  
          TSNOROT(ILMOS(K),JLMOS(K))=TSNOGAT(K)  
          RHOSROT(ILMOS(K),JLMOS(K))=RHOSGAT(K)  
          SNOROT (ILMOS(K),JLMOS(K))=SNOGAT (K)  
          GTROT  (ILMOS(K),JLMOS(K))=GTGAT  (K)                         
          WSNOROT(ILMOS(K),JLMOS(K))=WSNOGAT(K)  
          TCANROT(ILMOS(K),JLMOS(K))=TCANGAT(K)  
          RCANROT(ILMOS(K),JLMOS(K))=RCANGAT(K)  
          SCANROT(ILMOS(K),JLMOS(K))=SCANGAT(K)  
          GROROT (ILMOS(K),JLMOS(K))=GROGAT (K)  
          TACROT (ILMOS(K),JLMOS(K))=TACGAT (K)  
          QACROT (ILMOS(K),JLMOS(K))=QACGAT (K)  
          CMAIROT(ILMOS(K),JLMOS(K))=CMAIGAT(K)  
          REFROT (ILMOS(K),JLMOS(K))=REFGAT (K)                         
          BCSNROT(ILMOS(K),JLMOS(K))=BCSNGAT(K)                         
          EMISROT(ILMOS(K),JLMOS(K))=EMISGAT(K)                         
  100 CONTINUE
C
      DO 200 L=1,NBS                                                    
      DO 200 K=1,NML
          SALBROT(ILMOS(K),JLMOS(K),L)=SALBGAT(K,L)                     
          CSALROT(ILMOS(K),JLMOS(K),L)=CSALGAT(K,L)                     
  200 CONTINUE                                                          
C                                                                       
      DO 300 L=1,IG                                                     
      DO 300 K=1,NML                                                    
          TBARROT(ILMOS(K),JLMOS(K),L)=TBARGAT(K,L)
          THLQROT(ILMOS(K),JLMOS(K),L)=THLQGAT(K,L)
          THICROT(ILMOS(K),JLMOS(K),L)=THICGAT(K,L)
  300 CONTINUE                                                          
C
      DO 400 L=1,4                                                      
      DO 400 K=1,NML                                                    
          TSFSROT(ILMOS(K),JLMOS(K),L)=TSFSGAT(K,L)
  400 CONTINUE                                                          
C
      RETURN
      END
