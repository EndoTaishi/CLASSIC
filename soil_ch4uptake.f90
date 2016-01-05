subroutine soil_ch4uptake(Tbar,thpor,thliq,thice,bi,atm_CH4)

!           Canadian Terrestrial Ecosystem Model (CTEM)
!               Soil Methane Oxidation Subroutine
!
!History:

! >> J. Melton. Dec 22 2015 - Coded up based on C. Curry (2007) Modelling the
! the soil consumption of atmospheric methane at the global scale. Global
! Biogeo. Cycl. v. 21 GB4012 doi: 10.1029/2006GB002818.

use ctem_params,  only :

implicit none

! Arguments:
Tbar        ! Temperature of soil layers (K)
thpor       ! total porosity (cm^3 cm^-3)
bi          ! Clapp and Hornberger b-term (-)
atm_CH4         ! CH4 concentration at the soil surface (cm^-3) ? likely ppm CHECK!

! Local variables:
Tsoil       ! Temperature of soil layers (deg C)
real D_soil     ! Diffusion coefficient in soil (cm^2 s^-1)
G_T         ! Temperature factor used in determining D_soil
G_soil      ! Soil moisture factor used in determining D_soil
thpor_air    ! air-filled porosity (cm^3 cm^-3)
k_oxidr     ! first-order oxidation rate constant (s^-1)
k_o         ! base oxidation rate constant (s^-1)
r_T         ! temperature factor used in determination of rate constant
r_SM        ! soil moisture factor used in determination of rate constant

! Local parameters:
D_air = 0.196   ! Diffusion coefficient in air (cm^2 s^-1) @ STP
g_0 = 586.7 ! mg CH4 ppmv^-1 s d^-1 m^2

betaCH4 =

! ---------------------------------------------------------------------
! Begin

! Convert Tbar to Tsoil (from K to deg C)
Tsoil = Tbar - 273.16

! Find the diffusion coefficient in soil (D_soil)

G_T = 1.0 + 0.0055 * Tsoil


thpor_air = thpor - (Tliq + Tice) ! Check if already calc'd in CLASS.

! The bi  (Clapp and Hornberger b-term) is already calculated by CLASS as:
!bi = 15.9 * f_clay + 2.91 thus we use that value.

G_soil = thpor**(4/3) * (thpor_air / thpor)**(1.5 + 3 / bi)

D_soil = D_air * G_T * G_soil

! Determine the first-order oxidation rate constant (k_oxidr)
! first find the temperature term

if (Tsoil < 0.0 .and. Tsoil >= -10.0) then
    r_T = (0.1 * Tsoil + 1.0)**2
else if (Tsoil >= 0.0 .and. Tsoil < 43.3) then
    r_T = exp(0.0693 * Tsoil - 8.56E-7 * Tsoil**4)
else !all other temps (<-10 and >=43.3
    r_T = 0.
end if

! Next find the term based on soil moisture (suction)

! Find the soil water potential for the uppermost layer
! need the absolute value.
GRKINF(I,J)=ABS(GRKSATF(I,J)*(THLINF(I,J)/THPORF(I,J))
     1                        **(2.*BI(I,J)+3.))

      REAL GRKSAT(ILG,IG)   !Hydraulic conductivity of soil at
                            !saturation [m s-1] (Ksat) - need in kPa....
THLINF is calculated as the
      !maximum of finf(THPOR – THICE), THLIQ, and THLMIN, where finf
      !represents the fractional saturation of the soil behind the
      !wetting front, corresponding to a hydraulic conductivity of half
      !the saturated value (this correction is applied in order to
      !account for the fact that as infiltration occurs, a small amount
      !of air is usually trapped in the soil).

if ( GRKINF < 0.2) then !MPa!!
 r_SM = 1.0
else if (grkinf >= 0.2 .and. grkinf <= 100.) then !Mpa!!
 r_SM = (1. - (log10(grking) - log10(0.2)) / (log10(100) - log10(0.2)))**betaCH4
else ! grkinf > 100 MPa.
 r_SM = 0.
end if

k_oxidr = k_o * r_T * r_SM

! Find the flux correction for croplands

r_C =

! Find the flux correction for wetlands

r_W =

! Find the surface flux (J_ch4)

J_ch4 =  atm_CH4 * r_C * r_W * g_0 * (D_soil * k_oxidr)**0.5





! -------------
end subroutine soil_ch4uptake

