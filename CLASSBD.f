      BLOCK DATA 
C
C     Purpose: Assign values to parameters in CLASS common blocks. CLASS 
C     incorporates several kinds of parameters in its common blocks. 
C     Some are defined specifically for use in the CLASS code; some are 
C     also shared with the atmospheric model (if running in coupled 
C     mode).
C
C     * JAN 09/15 - E.CHAN/D.VERSEGHY. NEW VALUE FOR ALBRCK.
C     * MAR 13/09 - D.VERSEGHY. REPLACE SURFCON COMMON BLOCK WITH
C     *                         CLASSD2; NEW VARIABLE ANGMAX
C     * FEB 06/07 - D.VERSEGHY. NEW VALUE FOR ALVSO AND ALIRI.
C     * NOV 04/04 - D.VERSEGHY. NEW VALUES FOR ALVSI AND ALIRI.
C     * DEC 02/03 - D.VERSEGHY. HARMONIZE VALUES OF HCPICE,
C     *                         RHOICE AND SPHICE.
C     * JUL 18/03 - D.VERSEGHY. SPLIT OFF DATA STATEMENTS FROM
C     *                         OLD "CLASSD" ROUTINE TO COMPLY
C     *                         WITH FORTRAN 90 CONVENTION.
C     * AUG 06/02 - D.VERSEGHY. DEFINE PHYSICAL CONSTANTS PASSED
C     *                         THROUGH CLASS COMMON BLOCKS.  
C     *                         (BASED ON ROUTINE "HYDCON".)
C
C     * THE FOLLOWING COMMON BLOCK PARAMETERS ARE DEFINED WITHIN
C     * THE GCM.
C
      COMMON /PARAMS/ X1,    X2,    X3,    X4,   G,GAS,   X5,
     1                X6,    CPRES, GASV,  X7
      COMMON /PARAM1/ CPI,   X8,    CELZRO,X9,    X10,    X11
      COMMON /PARAM3/ X12,   X13,   X14,   X15,   SIGMA,  X16
      COMMON  /TIMES/ DELTIM,K1,    K2,    K3,    K4,     K5,
     1                K6,    K7,    K8,    K9,    K10,    K11

      COMMON /EPS/    A,B,EPS1,EPS2
      COMMON /HTCP/   T1S,T2S,AI,BI,AW,BW,SLP

C     * COMMON BLOCKS DEFINED SPECIFICALLY FOR USE IN CLASS.

      !
      !In routine CLASSBD, values are primarily assigned to the 
      !parameters that are specific to the CLASS code and are passed 
      !through it via common blocks CLASS1 through CLASS8. The table 
      !below lists the scalar parameters, their definitions, and their 
      !designated values with units.
      !
      !----------------------------------------------------------------|
      !Name   | Definition                    | Value      |   Units   |
      !----------------------------------------------------------------|
      !VKC    | Von Karman Constant           | 0.40       |   -       |
      !----------------------------------------------------------------|
      !CT     | Drag Coefficient for water    | 1.15*10^-3 |   -       |
      !----------------------------------------------------------------|
      !VMIN   | Minimum wind speed            | 0.1        |   m s-1   |
      !----------------------------------------------------------------|
      !TCW    | Thermal conductivity of water | 0.57       | W m-1 K-1 |
      !----------------------------------------------------------------|
      !TCICE  | Thermal conductivity of ice   | 2.24       | W m-1 K-1 |
      !----------------------------------------------------------------|
      !TCSAND | Thermal conductivity of sand  | 2.5        | W m-1 K-1 |
      !       | particles                     |            |           |
      !----------------------------------------------------------------|      
      !TCCLAY | Thermal conductivity of fine  | 2.5        | W m-1 K-1 |
      !       | mineral particles             |            |           |
      !----------------------------------------------------------------|
      !TCOM   | Thermal conductivity of       | 0.25       | W m-1 K-1 |
      !       | organic matter                |            |           |
      !----------------------------------------------------------------|
      !TCDRYS | Thermal conductivity of dry   | 0.275      | W m-1 K-1 |
      !       | mineral soil                  |            |           |      
      !----------------------------------------------------------------|
      !RHOSOL | Density of soil mineral matter| 2.65*10^3  | kg m-3    |
      !----------------------------------------------------------------|
      !RHOOM  | Density of soil organic matter| 1.30*10^3  | kg m-3    |
      !----------------------------------------------------------------|
      !HCPW   | Volumetric heat capacity of   | 4.187*10^6 | J m-3 K-1 |
      !       | water                         |            |           |
      !----------------------------------------------------------------|
      !HCPICE | Volumetric heat capacity of   | 1.9257*10^6| J m-3 K-1 |
      !       | ice                           |            |           |
      !----------------------------------------------------------------|
      !HCPSOL | Volumetric heat capacity of   | 2.25*10^6  | J m-3 K-1 | 
      !       | mineral matter                |            |           |
      !----------------------------------------------------------------|
      !HCPOM  | Volumetric heat capacity of   | 2.50*10^6  | J m-3 K-1 |
      !       | organic matter                |            |           |
      !----------------------------------------------------------------|
      !HCPSND | Volumetric heat capacity of   | 2.13*10^6  | J m-3 K-1 |
      !       | sand particles                |            |           |
      !----------------------------------------------------------------|
      !HCPCLY | Volumetric heat capacity of   | 2.38*10^6  | J m-3 K-1 |
      !       | fine mineral particles        |            |           |
      !----------------------------------------------------------------|
      !SPHW   | Specific heat of water        | 4.186*10^3 | J kg-1 K-1|
      !----------------------------------------------------------------|
      !SPHICE | Specific heat of ice          | 2.10*10^3  | J kg-1 K-1|
      !----------------------------------------------------------------|
      !SPHVEG | Specific heat of vegetation   | 2.70*10^3  | J kg-1 K-1|
      !       | matter                        |            |           |
      !----------------------------------------------------------------|
      !RHOW   | Density of water              | 1.0*10^3   | kg m-3    |
      !----------------------------------------------------------------|
      !RHOICE | Density of ice                | 0.917*10^3 | kg m-3    |
      !----------------------------------------------------------------|
      !TCGLAC | Thermal conductivity of ice   | 2.24       | W m-1 K-1 |
      !       | sheets                        |            |           |
      !----------------------------------------------------------------|
      !CLHMLT | Latent heat of freezing of    | 0.334*10^6 | J kg-1    |
      !       | water                         |            |           |
      !----------------------------------------------------------------|
      !CLHVAP | Latent heat of vaporization   | 2.501*10^6 | J kg-1    |
      !       | of water                      |            |           |
      !----------------------------------------------------------------|
      !ZOLNG  | Natural log of roughness      | -4.605     | -         |
      !       | length of soil                |            |           |
      !----------------------------------------------------------------|
      !ZOLNS  | Natural log of roughness      | -6.908     | -         |
      !       | length of snow                |            |           |
      !----------------------------------------------------------------|
      !ZOLNI  | Natural log of roughness      | -6.215     | -         |
      !       | length of ice                 |            |           |
      !----------------------------------------------------------------|
      !ZORATG | Ratio of soil roughness length| 3.0        | -         |
      !       | for momentum to roughness     |            |           |
      !       | length for heat               |            |           |
      !----------------------------------------------------------------|
      !ALVSI  | Visible albedo of ice         | 0.95       | -         |
      !----------------------------------------------------------------|
      !ALIRI  | Near-infrared albedo of ice   | 0.73       | -         |
      !----------------------------------------------------------------|
      !ALVSO  | Visible albedo of organic     | 0.05       | -         |
      !       | matter                        |            |           |
      !----------------------------------------------------------------|
      !ALIRO  | Near-infrared albedo of       | 0.30       | -         |
      !       | organic matter                |            |           |
      !----------------------------------------------------------------|
      !ALBRCK | Albedo of rock                | 0.27       | -         |
      !----------------------------------------------------------------|
      !
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
     1                RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
      COMMON /CLASS5/ THPORG(3),THRORG(3),THMORG(3),BORG(3),
     1                PSISORG(3),GRKSORG(3)
      COMMON /CLASS6/ PI,GROWYR(18,4,2),ZOLNG,ZOLNS,ZOLNI,ZORAT(4),
     1                ZORATG
      COMMON /CLASS7/ CANEXT(4),XLEAF(4)
      COMMON /CLASS8/ ALVSI,ALIRI,ALVSO,ALIRO,ALBRCK
      COMMON /PHYCON/ DELTA,CGRAV,CKARM,CPD
      COMMON /CLASSD2/ AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX
      COMMON /ESTWI/  RW1,RW2,RW3,RI1,RI2,RI3
C
      DATA      VKC,        CT,         VMIN
     1       /  0.40,       1.15E-3,    0.1     /

      DATA      TCW,        TCICE,      TCSAND,     TCCLAY,     TCOM
     1       /  0.57,       2.24,       2.5,        2.5,        0.25   /
C     1       /  0.57,       2.24,       8.0,        2.5,        0.25   /

      DATA      TCDRYS,     RHOSOL,     RHOOM
     1       /  0.275,      2.65E3,     1.30E3  /    

      DATA      HCPW,       HCPICE,     HCPSOL,     HCPOM
     1       /  4.187E6,    1.9257E6,   2.25E6,     2.50E6 /

      DATA      HCPSND,     HCPCLY,     SPHW,       SPHICE,     SPHVEG
     1       /  2.13E6,     2.38E6,     4.186E3,    2.10E3,     2.70E3 /

      DATA      RHOW,       RHOICE,     TCGLAC,     CLHMLT,     CLHVAP
     1       /  1.0E3,      0.917E3,    2.24,       0.334E6,    2.501E6/

      DATA      ZOLNG,      ZOLNS,      ZOLNI,      ZORATG
     1       /  -4.605,     -6.908,     -6.215,     3.0         /

      DATA      ALVSI,      ALIRI,      ALVSO,      ALIRO,      ALBRCK
     1       /  0.95,       0.73,       0.05,       0.30,       0.16  /

      !Values are also assigned to several non-scalar parameters, as 
      !follows:
      !
      !1)   The crop growth descriptor array GROWYR (see the 
      !     documentation for subroutine APREP);
      !2)   Three parameters for the four main vegetation categories 
      !     recognized by CLASS (needleleaf trees, broadleaf trees, 
      !     crops and grass): ZORAT, the ratio of the roughness length 
      !     for momentum to the roughness length for heat (currently set 
      !     to 1); CANEXT, an attenuation coefficient used in 
      !     calculating the sky view factor for vegetation canopies 
      !     (variable c in the documentation for subroutine CANALB); and 
      !     XLEAF, a leaf dimension factor used in calculating the leaf 
      !     boundary resistance (variable Cl in the documentation for 
      !     subroutine APREP);
      !
      !3)   Six hydraulic parameters associated with the three basic 
      !     types of organic soils (fibric, hemic and sapric): THPORG, 
      !     THRORG, THMORG, BORG, PSISORG and GRKSORG (see the 
      !     documentation for subroutine CLASSB). The table below lists 
      !     their values, derived from the work of Letts et al. (2000), 
      !     alongside the symbols used in the CLASSB documentation.
      !
      !----------------------------------------------------------------|
      !Name    | Symbol    | Fibric peat  | Hemic peat  | Sapric peat  |
      !----------------------------------------------------------------|
      !THPORG  | theta_p   | 0.93         | 0.88        | 0.83         |
      !----------------------------------------------------------------|
      !BORG    | b         | 2.7          | 6.1         | 12.0         |
      !----------------------------------------------------------------|
      !GRKSORG | K_sat     | 2.8*10^-4    | 2.0*10^-6   | 1.0*10^-7    |
      !----------------------------------------------------------------|
      !PSISORG | psi_sat   | 0.0103       | 0.0102      | 0.0101       |
      !----------------------------------------------------------------|
      !THMORG  | theta_min | 0.04         | 0.15        | 0.22         |
      !----------------------------------------------------------------|
      !THRORG  | theta_ret | 0.275        | 0.62        | 0.705        |
      !----------------------------------------------------------------|
      !
      DATA GROWYR   /213.,213.,213.,213.,213.,213.,0.,0.,0.,
     1               0.,0.,0., 75.,106.,136.,167.,167.,167.,
     2               273.,273.,273.,273.,273.,273.,0.,0.,0.,
     3               0.,0.,0.,135.,166.,196.,196.,196.,196.,
     4               121.,121.,121.,121.,121.,121.,0.,0.,0.,
     5               0.,0.,0.,275.,244.,214.,214.,214.,214.,
     6               151.,151.,151.,151.,151.,151.,0.,0.,0.,
     7               0.,0.,0.,305.,274.,244.,244.,244.,244.,
     8               213.,213.,213.,213.,213.,213.,0.,0.,0.,
     9               0.,0., 75.,106.,136.,167.,167.,167.,167.,
     A               273.,273.,273.,273.,273.,273.,0.,0.,0.,
     B               0.,0.,135.,166.,196.,196.,196.,196.,196.,
     C               121.,121.,121.,121.,121.,121.,0.,0.,0.,
     D               0.,0.,275.,244.,214.,214.,214.,214.,214.,
     E               151.,151.,151.,151.,151.,151.,0.,0.,0.,
     F               0.,0.,305.,274.,244.,244.,244.,244.,244. /
C
      DATA ZORAT    /1.0,1.0,1.0,1.0/
      DATA CANEXT   /-0.5,-1.5,-0.8,-0.8/
      DATA XLEAF    /0.0247,0.0204,0.0456,0.0456/
      DATA THPORG   /0.93,0.88,0.83/
      DATA THRORG   /0.275,0.620,0.705/
      DATA THMORG   /0.04,0.15,0.22/
      DATA BORG     /2.7,6.1,12.0/
      DATA PSISORG  /0.0103,0.0102,0.0101/
      DATA GRKSORG  /2.8E-4,2.0E-6,1.0E-7/
      !
      !Finally, if CLASS is being run offline, values must be assigned 
      !to the parameters listed in the first table which would normally 
      !be assigned in the GCM:
      !
      !----------------------------------------------------------------|
      !GCM Name | Definition                  | Value       | Units    |
      !----------------------------------------------------------------|
      !DELTIM   | Time Step                   |Varies by run| s        |
      !----------------------------------------------------------------|
      !CELZRO   | Freezing point of water     | 273.16      | K        |    
      !----------------------------------------------------------------|
      !GAS      | Gas constant                | 287.04      |J kg-1 K-1|
      !----------------------------------------------------------------|
      !GASV     |Gas constant for water vapour| 461.50      |J kg-1 K-1| 
      !----------------------------------------------------------------|
      !G        | Acceleration due to gravity | 9.80616     | m s-1    |   
      !----------------------------------------------------------------|
      !SIGMA    | Stefan-Boltzmann constant   |5.66796*10^-8| W m-2 K-4|
      !----------------------------------------------------------------|
      !CPRES    | Specific heat of air        | 1.00464*10^3|J kg-1 K-1|
      !----------------------------------------------------------------|
      !CPI      | Pi                          | 3.14159265  | -        |
      !----------------------------------------------------------------|
      !
C
C     * ASSIGN VALUES NORMALLY SPECIFIED WITHIN THE GCM.
C
      DATA      G,          GAS,        CPRES,      GASV
     1/         9.80616,    287.04,     1.00464E3,  461.50   /

      DATA      CELZRO,     SIGMA,      DELTIM
C     1/         273.16,     5.66796E-8, 900.0   /
     1/         273.16,     5.66796E-8, 1800.0   /

      DATA      CPI 
     1       /  3.1415926535898    /
      !
      !The parameters in common blocks PHYCON and CLASSD2 that do not 
      !have corresponding values assigned in the GCM are assigned their 
      !RPN values, and the remaining parameters in the GCM common blocks 
      !PARAMS, PARAM1, PARAM3 and TIMES are assigned dummy values.
      !
C
C     * ADDITIONAL VALUES FOR RPN AND GCM COMMON BLOCKS.
C
      DATA      DELTA,      AS,         ASX,        ANGMAX
     1       /  0.608,      12.0,       4.7,        0.85   /

      DATA      CI,         BS,         BETA,       FACTN,      HMIN 
     1       /  40.0,       1.0,        1.0,        1.2,        40.   /

      DATA      A,          B,          EPS1,       EPS2
     1       /  21.656,     5418.0,     0.622,      0.378  /

      DATA      T1S,        T2S
     1       /  273.16,     233.16   /

      DATA      RW1,        RW2,        RW3
     1       /  53.67957,  -6743.769,  -4.8451   /

      DATA      RI1,        RI2,        RI3
     1       /  23.33086,  -6111.72784, 0.15215   /

      DATA      X1,  X2,  X3,  X4,  X5,  X6,  X7,  X8
     1       /  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  /

      DATA      X9,  X10, X11, X12, X13, X14, X15, X16
     1       /  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  /

      DATA      K1,  K2,  K3,  K4,  K5,  K6,  K7,  K8
     1       /  0,   0,   0,   0,   0,   0,   0,   0    /

      DATA      K9,  K10, K11
     1       /  0,   0,   0    /
C
      END
