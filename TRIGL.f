      SUBROUTINE TRIGL(ILATH,SR,WR,CR,RADR,WOSQ)
C 
C     * JUL 14/92 - E. CHAN (ADD REAL*8 DECLARATIONS)
C     * JAN 19/78 - J.D.HENDERSON 
C     * THE ARGUMENT LIST IS THE SAME AS FOR GAUSSG.
C     * GAUSSG FILLS ONLY THE N HEM ORDERED N TO S. 
C     * THIS ROUTINE MAKES THE ARRAYS GLOBAL AND ORDERED FROM S TO N. 
C     *      SR=SIN(LAT),  CR=COS(LAT),  RADR=LATITUDE IN RADIANS.
C     *      WR = GAUSSIAN WEIGHTS,  WOSQ = WR/(SR**2). 

      use ctem_params, only: lat
C 
      REAL*8 SR(lat),WR(lat),CR(lat)
      REAL*8 RADR(lat),WOSQ(lat)
C-------------------------------------------------------------------- 
C     * CR,WR,WOSQ ARE SYMMETRIC ABOUT THE EQUATOR. 
C     * SR AND RAD ARE ANTISYMMETRIC. 
C 
      PIH=3.14159265/2. 
      ILAT=ILATH*2
      DO 150 J=1,ILATH
      K=ILAT+1-J
      CR(K)=CR(J) 
      WR(K)=WR(J) 
      WOSQ(K)=WOSQ(J) 
      SR(K)=SR(J) 
      SR(J)=-SR(J)
      RADR(K)=PIH-RADR(J) 
      RADR(J)=-RADR(K)
  150 CONTINUE
C 
      RETURN
      END 
