subroutine soil_ch4uptake(IL1,IL2,TBAR,THP,BI,THLQ,THIC, &
     &                     PSIS,GRAV,FCAN,obswetf,wetfdyn, &
     &                     wetfracgrd,SAND,RHOW,RHOICE,CH4_soills)
!           Canadian Terrestrial Ecosystem Model (CTEM)
!               Soil Methane Oxidation Subroutine
!
!History:

! >> J. Melton. Dec 22 2015 - Coded up based on C. Curry (2007) Modelling the
! the soil consumption of atmospheric methane at the global scale. Global
! Biogeo. Cycl. v. 21 GB4012 doi: 10.1029/2006GB002818.

use ctem_params,  only : ilg,ignd,icp1,nlat,wtCH4

implicit none

! Arguments:
real, dimension(ilg,ignd), intent(in) :: TBAR       ! Temperature of soil layers (K)
real, dimension(ilg,ignd), intent(in) :: THP        ! Total porosity ($cm^3 cm^{-3}$)
real, dimension(ilg,ignd), intent(in) :: BI         ! Clapp and Hornberger b-term (-)
real, dimension(ilg,ignd), intent(in) :: THLQ       ! Fractional water content (-)
real, dimension(ilg,ignd), intent(in) :: THIC       ! Fractional ice content (-)
real, dimension(ilg,ignd), intent(in) :: PSIS       ! Soil moisture suction at saturation (m)
real, dimension(ilg,ignd), intent(in) :: SAND       ! Sand content (%)
real, dimension(ilg,icp1), intent(in) :: FCAN       ! Fractional coverage of vegetation (-)
real, dimension(nlat), intent(in) :: wetfracgrd     ! Prescribed fraction of wetlands in a grid cell
real, dimension(ilg), intent(in) :: wetfdyn         ! Dynamic gridcell wetland fraction determined using slope and soil moisture
real, dimension(ilg), intent(in) :: atm_CH4         ! CH4 concentration at the soil surface ($cm^-3$) ? likely ppm CHECK!
real, intent(in) :: GRAV                            ! Acceleration due to gravity ($m s^{-1}$)
logical, intent(in) :: obswetf                      ! Switch, if true then use the prescribed wetland cover
real, intent(in) :: RHOW                            ! Density of water ($kg m^{-3}$)
real, intent(in) :: RHOICE                          ! Density of ice ($kg m^{-3}$)

real, dimension(ilg), intent(out) :: CH4_soills          ! Methane uptake into the soil column ($mg CH_4 m^{-2} s^{-1}$)

! Local variables:
real :: Tsoil                                       ! Temperature of soil layers ($/circ$C)
real :: D_soil                                      ! Diffusivity of CH4 in soil ($cm^2 s^{-1}$)
real :: G_T                                         ! Temperature factor used in determining D_soil (-)
real :: G_soil                                      ! Soil moisture factor used in determining D_soil (-)
real :: THP_air                                     ! Air-filled porosity ($cm^3 cm^{-3}$)
real :: k_oxidr                                     ! First-order oxidation rate constant ($s^-1$)
real :: r_T                                         ! Temperature factor used in determination of rate constant (-)
real :: r_SM                                        ! Soil moisture factor used in determination of rate constant (-)
integer :: i,j,layer                                ! Counters
integer :: isand                                    ! Check if bedrock (value of -3)

! Local parameters:
real, parameter :: D_air = 0.196        ! Diffusivity of CH4 in air (cm^2 s^-1) @ STP
real, parameter :: g_0 = 586.7 / 86400. ! Scaling factor takes CH4_soills to mg CH4 m^-2 s^-1 (units: mg $CH_4 ppmv^{-1} s s^{-1} m^{-2} cm{-1}$)
real, parameter :: betaCH4 = 0.8        ! Constant derived in Curry (2007) from comparison against measurements (-)
real, parameter :: k_o = 5.03E-5        ! Base oxidation rate derived in Curry (2007) from comparison against measurements ($s^{-1}$)

! ---------------------------------------------------------------------
! Begin

! Find the flux correction due to wetlands

if (obswetf) then  ! Use the prescribed wetland fractions
    r_W = 1.0 - wetfracgrd
else ! use the dynamically determined wetland area
    r_W = 1.0 - wetfdyn
end if

! The soil oxidation methane sink is assumed to only operate in the first model
! soil layer, thus we only consider that layer here.
layer = 1

do 10 i = il1, il2

    isand = nint(SAND(i,layer))

    if (isand <= -1) goto 10 ! not soil so move on.

    ! Convert TBAR to Tsoil (from K to deg C)
    Tsoil = TBAR(i,layer) - 273.16

    ! Find the diffusion coefficient in soil (D_soil)

    ! First the temperature factor, G_T:
    G_T = 1.0 + 0.0055 * Tsoil

    ! Find the air filled porosity, THP_air:
    THP_air = THP(i,layer) - (THLQ(i,layer) + THIC(i,layer)*RHOICE/RHOW)
    ! Note: THP_air can fall to < 0 after snow melt
    if (THP_air  < 0.) then
        write(*,*)'In soil_ch4uptake THP_air became negative!'
        call exit()
        !THP_air = max(0.,THP_air)
    end if

    ! The BI  (Clapp and Hornberger b-term) is already calculated by CLASS as:
    !BI = 15.9 * f_clay + 2.91, thus we use that value.

    ! G_soil is the influence of the soil texture, moisture, and porosity:
    G_soil = THP(i,layer)**(4./3.) * (THP_air / THP(i,layer))**(1.5 + 3. / BI(i,layer))

    ! The diffusion coefficient of CH4 in soil is then:
    D_soil = D_air * G_T * G_soil

    ! Determine the first-order oxidation rate constant (k_oxidr)

    ! First find the temperature term, r_T (FLAG note that Charles' original code does not have the high temp limit!)

    if (Tsoil < 0.0 .and. Tsoil >= -10.0) then
        r_T = (0.1 * Tsoil + 1.0)**2
    else if (Tsoil >= 0.0 .and. Tsoil < 43.3) then
        r_T = exp(0.0693 * Tsoil - 8.56E-7 * Tsoil**4)
    else !all other temps (<-10 and >=43.3)
        r_T = 0.
    end if

    ! Next find the term based on soil moisture (suction)

    ! Find the soil water potential for the uppermost layer
    ! need the absolute value.
    psi = abs(PSIS(i,layer) * (THLQ(i,layer)/THP(i,layer))**(-BI(i,layer)))

    ! Convert units from m to kPa
    psi = psi * GRAV

    if ( psi < 200.) then !0.2 MPa in paper (NOTE: In Charles's code this is <=, but is < in paper)
        r_SM = 1.0
    else if (psi >= 200. .and. psi <= 100000.) then !0.2 and 100 Mpa in paper
        r_SM = (1. - (log10(psi) - log10(0.2)) / (log10(100) - log10(0.2)))**betaCH4
    else ! psi > 100 MPa.
        r_SM = 0.
    end if

    k_oxidr = k_o * r_T * r_SM

    ! Find the flux correction for croplands

    r_C = 1.0 - (0.75 * FCAN(i, 3))

    ! Find the surface flux (CH4_soills) for each tile, then for each gridcell

    CH4_soills(i) =  atm_CH4 * r_C * r_W * g_0 * sqrt(D_soil * k_oxidr)

    ! Convert from mg CH4 m^-2 s^-1 to umol CH4 m^-2 s^-1

    CH4_soills(i) = CH4_soills(i) * 1.E3 / wtCH4


10 continue

end subroutine soil_ch4uptake

