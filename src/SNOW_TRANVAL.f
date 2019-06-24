!>\file
C!Computes the direct and diffuse snow transmission using lookup table and current snow conditions.
C!@author J. Cole
!!
C!! This subroutine compute the direct and diffuse snow transmission using a
!! lookup table and information about the current snow pack state.
!! Transmission are computed for each solar radiation wavelength intervals
!! so a total of 8 albedos will be returned.  These transmissions can then be
!! used to compute the total snow transmission based on the by weighting
!! the results by the direct beam fraction of the incident solar radiation.
C!
      SUBROUTINE SNOW_TRANVAL(trandif, ! OUTPUT
     1                        trandir,
     2                        smu,     ! INPUT
     3                        salb,
     4                        bc_conc,
     5                        snow_reff,
     6                        swe,
     7                        c_ind,
     8                        il1,
     9                        il2,
     1                        ilg,
     2                        nbnd     )
!
!     * Feb 10/2015 - J.Cole. New version for gcm18:
!                             - NBC increased from 12 to 20.
!                               Therefor LBC_CONC data statement
!                               changed accordingly.
!     * JAN 24/2013 - J.COLE. Previous version for gcm17:
!                    - COMPUTES THE DIRECT AND DIFFUSE SNOW TRANSMISSION
!                      USING LOOKUP TABLE AND CURRENT SNOW CONDITIONS.
!
      use classic_params, only : nsmu,nsalb,nbc,nreff,nswe,nbnd_lut,
     1                        trandif_lut,trandir_lut
      IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! THIS SUBROUTINE COMPUTE THE DIRECT AND DIFFUSE SNOW TRANSMISSION USING A
! LOOKUP TABLE AND INFORMATION ABOUT THE CURRENT SNOW PACK STATE.  
! TRANSMISSION ARE COMPUTED FOR EACH SOLAR RADIATION WAVELENGTH INTERVALS
! SO A TOTAL OF 8 ALBEDOS WILL BE RETURNED.  THESE TRANSMISSIONS CAN THEN BE
! USED TO COMPUTE THE TOTAL SNOW TRAMISSION BY WEIGHTING
! THE RESULTS BY THE DIRECT BEAM FRACTION OF THE INCIDENT SOLAR RADIATION.
!
! INPUTS
! SMU:       COSINE OF THE SOLAR ZENITH ANGLE [UNITLESS]
! SALB :     ALBEDO OF THE UNDERLYING SURFACE [UNITLESS]
! BC_CONC:   CONCENTRATION OF BLACK CARBON IN THE SNOW PACK [NG (BC)/KG (SNOW)]
! SNOW_REFF: EFFECTIVE RADIUS OF THE SNOW GRAIN [MICRONS]
! SWE:       SNOW WATER EQUIVALENT (SNOWPACK DENSITY*SNOW PACK DEPTH) [KG/M^2]
! C_IND:     INDICATOR THAT A CALCULATION SHOULD BE PERFORMED FOR THIS POINT
!            1-YES, 0-NO
! IL1:       STARTING POINT FOR TRANSMISSION CALCULATIONS
! IL2:       ENDING POINT FOR TRANSMISSION CALCULATIONS
! ILG:       NUMBER OF POINTS FOR WHICH TO COMPUTE TRANSMISSIONS
! NBND:      NUMBER OF WAVELENGTH INTERVALS FOR WHICH TO COMPUTE THE TRANSMISSIONS
!
! OUTPUTS
! TRANDIF: DIFFUSE SNOW TRANSMISSION (AKA WHITE SKY TRANSMISSION)
! TRANDIR: DIRECT BEAM SNOW TRANSMISSION (AKA BLACK SKY TRANSMISSION)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!
! INPUT
!
      REAL, INTENT(IN), DIMENSION(ILG) ::
     1 smu,
     2 bc_conc,
     3 snow_reff,
     4 swe

      REAL, INTENT(IN), DIMENSION(ILG,NBND) ::
     1 salb
 
      INTEGER, INTENT(IN), DIMENSION(ILG) ::
     1 c_ind

      INTEGER, INTENT(IN) ::
     1 il1,
     2 il2,
     3 ilg,
     4 nbnd

!
! OUTPUT
!
      REAL, INTENT(OUT), DIMENSION(ILG,NBND) ::
     1 trandif,
     2 trandir

!
! LOCAL
!
      REAL, DIMENSION(ILG,2) ::
     1 wsmu,
     2 wbc,
     3 wreff,
     4 wswe

      REAL :: 
     1 wsalb(2)

      REAL :: 
     1 wtt

      INTEGER, DIMENSION(ILG) :: 
     1 ismu,
     2 ibc,
     3 ireff,
     4 iswe

      INTEGER ::
     1 ib,
     2 i,
     3 isalb

      INTEGER ::
     1 iismu,
     2 iisalb,
     3 iibc,
     4 iireff,
     5 iiswe

      INTEGER ::
     1 mvidx
     
!
! CONSTANTS
!
!      INTEGER, PARAMETER ::
!     1 nsmu     = 10,
!     2 nsalb    = 11,
!     3 nbc      = 20,
!     4 nreff    = 10,
!     5 nswe     = 11,
!     6 nbnd_lut = 4

      REAL, PARAMETER :: ! STATE VALUES FOR LUT
     1 LSALB(NSALB)    = 
     2                  (/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/),
     3 LSMU(NSMU)        = 
     4                      (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/),
     5 LSNOW_REFF(NREFF) = 
     6   (/50.0,75.0,100.0,150.0,200.0,275.0,375.0,500.0,700.0,1000.0/),
     7 LSWE(NSWE)        = 
     8      (/0.1,0.25,0.65,1.7,4.4,12.0,30.0,75.0,200.0,500.0,5000.0/),
     9 LBC_CONC(NBC)     = (/0.0,1.0,
     1                       5.0,10.0,
     2                       50.0,100.0,
     3                       500.0,1000.0,
     4                       5000.0,10000.0,
     5                       50000.0,100000.0, 
     6                       250000.0,500000.0,750000.0,1000000.0,
     7                       2500000.0,5000000.0,7500000.0,10000000.0/)

!      REAL, DIMENSION(NBC,NSWE,NREFF,NSMU,NSALB,NBND_LUT) ::
!     1 trandif_lut,
!     2 trandir_lut

!      INTEGER ::
!     1 snow_tran_lut_init
      
!       COMMON /SNOWTRANLUT/ trandif_lut,trandir_lut,snow_tran_lut_init

! ABORT IF THE LUT HAS NOT BEEN READ IN
!      IF (snow_tran_lut_init .NE. 1) THEN
!         WRITE(6,*) 'SNOW TRANSMISSION LUT HAS NOT BEEN INITIALIZED',
!     1              snow_tran_lut_init
!         CALL XIT('SNOW_TRANVAL',-1)
!      END IF

! ABORT IF THE NUMBER OF BANDS IN THE LUT DOES NOT MATCH SIZE PASSED IN
      IF (nbnd_lut .NE. nbnd) THEN
         WRITE(6,*) 'MISMATCH IN NUMBER OF WAVELENGTH INTERVALS'
         CALL XIT('SNOW_TRANVAL',-2)
      END IF

! COMPUTE THE TRANSMISSIONS USING LINEAR INTERPOLATION

! COMPUTE THE INTERPOLATION WEIGHTS AND POINTS ONCE AND REUSE FOR
! TRANSMISSION INTERPOLATION FOR EACH BAND.
! HAVE A CHECK TO SET THE WEIGHTS DEPENDING IF THE INPUT IS
! OUTSIDE OR INSIDE THE LOOKUP TABLE RANGE
      
      DO i = il1, il2
         IF (c_ind(i) .EQ. 1) THEN
            ismu(i)  = mvidx(LSMU,       nsmu,  smu(i))
            ibc(i)   = mvidx(LBC_CONC,   nbc,   bc_conc(i))
            ireff(i) = mvidx(LSNOW_REFF, nreff, snow_reff(i))
            iswe(i)  = mvidx(LSWE,       nswe,  swe(i))

            IF (smu(i) .LE. LSMU(1)) THEN
               wsmu(i,2) = 0.0
               wsmu(i,1) = 1.0-wsmu(i,2)
            ELSEIF (smu(i) .GT. LSMU(NSMU)) THEN
               wsmu(i,2) = 1.0
               wsmu(i,1) = 1.0-wsmu(i,2)
            ELSE
               wsmu(i,2) = (smu(i)-LSMU(ismu(i)))
     1                   / (LSMU(ismu(i)+1)-LSMU(ismu(i)))
               wsmu(i,1) = 1.0-wsmu(i,2)
            END IF
            
            IF (bc_conc(i) .LE. LBC_CONC(1)) THEN
               wbc(i,2) = 0.0
               wbc(i,1) = 1.0-wbc(i,2)
            ELSEIF (bc_conc(i) .GT. LBC_CONC(NBC)) THEN
               wbc(i,2) = 1.0
               wbc(i,1) = 1.0-wbc(i,2)            
            ELSE
               wbc(i,2) = (bc_conc(i)-LBC_CONC(ibc(i)))
     1                  / (LBC_CONC(ibc(i)+1)-LBC_CONC(ibc(i)))
               wbc(i,1) = 1.0-wbc(i,2)
            END IF
            
            IF (snow_reff(i) .LE. LSNOW_REFF(1)) THEN
               wreff(i,2) = 0.0
               wreff(i,1) = 1.0-wreff(i,2)
            ELSEIF (snow_reff(i) .GT. LSNOW_REFF(NREFF)) THEN
               wreff(i,2) = 1.0
               wreff(i,1) = 1.0-wreff(i,2)
            ELSE
               wreff(i,2) = (snow_reff(i)-LSNOW_REFF(ireff(i)))
     1                    / (LSNOW_REFF(ireff(i)+1)
     2                                   -LSNOW_REFF(ireff(i)))
               wreff(i,1) = 1.0-wreff(i,2)
            END IF
            
            IF (swe(i) .LE. LSWE(1)) THEN
               wswe(i,2) = 0.0
               wswe(i,1) = 1.0-wswe(i,2)
            ELSEIF (swe(i) .GT. LSWE(NSWE)) THEN
               wswe(i,2) = 1.0
               wswe(i,1) = 1.0-wswe(i,2)
            ELSE
               wswe(i,2) = (swe(i)-LSWE(iswe(i)))
     1                   / (LSWE(iswe(i)+1)-LSWE(iswe(i)))
               wswe(i,1) = 1.0-wswe(i,2)
            END IF
         END IF
      END DO ! i

      DO ib = 1, nbnd
         DO i = il1, il2
            IF (c_ind(i) .EQ. 1) THEN

               isalb = mvidx(LSALB,    nsalb, salb(i,ib))
               
               IF (salb(i,ib) .LE. LSALB(1)) THEN
                  wsalb(2) = 0.0
                  wsalb(1) = 1.0-wsalb(2)
               ELSEIF (salb(i,ib) .GT. LSALB(NSALB)) THEN
                  wsalb(2) = 1.0
                  wsalb(1) = 1.0-wsalb(2)
               ELSE
                  wsalb(2) = (salb(i,ib)-LSALB(isalb))
     1                     / (LSALB(isalb+1)-LSALB(isalb))
                  wsalb(1) = 1.0-wsalb(2)
               END IF
               
               trandir(i,ib) = 0.0
               trandif(i,ib) = 0.0

               DO iisalb = isalb,isalb+1
                  DO iismu = ismu(i),ismu(i)+1
                     DO iireff = ireff(i),ireff(i)+1
                        DO iiswe = iswe(i), iswe(i)+1
                           DO iibc = ibc(i), ibc(i)+1
                              
                              wtt = wsmu(i,iismu-ismu(i)+1)
     +                            * wreff(i,iireff-ireff(i)+1)
     +                            * wswe(i,iiswe-iswe(i)+1)
     +                            * wbc(i,iibc-ibc(i)+1)
     +                            * wsalb(iisalb-isalb+1)                              

                              trandif(i,ib) = trandif(i,ib) + wtt
     +                   *trandif_lut(iibc,iiswe,iireff,iismu,iisalb,ib)
                              trandir(i,ib) = trandir(i,ib) + wtt
     +                   *trandir_lut(iibc,iiswe,iireff,iismu,iisalb,ib)

                           END DO ! iibc
                        END DO  ! iiswe
                     END DO     ! iireff
                  END DO        ! iismu
               END DO           ! iisalb

               IF(trandif(i,ib) .GT. 1.3 .OR. 
     1                                     trandif(i,ib) .LT. 0.0) THEN
                  WRITE(6,*) 'Bad trandif ',i,ib,smu(i),bc_conc(i),
     1                     snow_reff(i),swe(i),salb(i,ib),trandif(i,ib)
                  WRITE(6,*) i,ib,ismu(i),ibc(i),ireff(i),iswe(i),isalb
                  CALL XIT('SNOW_TRANVAL',-3)
               END IF
               IF(trandir(i,ib) .GT. 1.3 .OR. 
     1                                     trandir(i,ib) .LT. 0.0) THEN
                  WRITE(6,*) 'Bad trandir ',i,ib,smu(i),bc_conc(i),
     1                      snow_reff(i),swe(i),salb(i,ib),trandir(i,ib)
                  WRITE(6,*) i,ib,ismu(i),ibc(i),ireff(i),iswe(i),isalb
                  CALL XIT('SNOW_TRANVAL',-3)
               END IF
            ELSE
               trandif(i,ib) = -999.0
               trandir(i,ib) = -999.0
            END IF
         END DO                 ! i
      END DO                    ! ib
      
      RETURN
      !>\file
      END                                                                         SNOW_TRANVAL.288
