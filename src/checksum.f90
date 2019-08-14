!> \file
!> Runs checksum calculations across several model variables
!!
module checksum

  implicit none

  public :: checksumCalc
  public :: bitcount, bitcount_int

contains
  !> \ingroup checksum_checksumCalc
  !! @{
  !! This subroutine takes the lonIndex and latIndex of a cell, accesses many
  !! attributes of the cell after the run, and creates a checksum from those attributes.
  !! This checksum is written to a .csv file in the output directory, which is then
  !! compared against the checksum from a previous run.
  !! @author M. Fortier
  !!
  subroutine checksumCalc(lonIndex, latIndex)

    use ctem_statevars, only : c_switch, vrot
    use class_statevars, only : class_rot
    use classic_params, only : icc, nmos, ignd, icp1, modelpft, iccp2, TFREZ

    implicit none

    ! arguments
    integer, intent(in) :: lonIndex       !< Index of the longitude of this cell
    integer, intent(in) :: latIndex       !< Index of the latitude of this cell
    ! pointers:
    real, pointer, dimension(:,:,:) :: FCANROT
    real, pointer, dimension(:,:)   :: FAREROT
    real, pointer, dimension(:,:,:) :: TBARROT
    real, pointer, dimension(:,:,:) :: THLQROT
    real, pointer, dimension(:,:,:) :: THICROT
    real, pointer, dimension(:,:)   :: TCANROT
    real, pointer, dimension(:,:)   :: TSNOROT
    real, pointer, dimension(:,:)   :: TPNDROT
    real, pointer, dimension(:,:)   :: ZPNDROT
    real, pointer, dimension(:,:)   :: RCANROT
    real, pointer, dimension(:,:)   :: SCANROT
    real, pointer, dimension(:,:)   :: SNOROT
    real, pointer, dimension(:,:)   :: ALBSROT
    real, pointer, dimension(:,:)   :: RHOSROT
    real, pointer, dimension(:,:)   :: GROROT

    logical, pointer :: ctem_on
    logical, pointer :: PFTCompetition
    logical, pointer :: lnduseon
    real, pointer, dimension(:,:,:) :: fcancmxrow
    real, pointer, dimension(:,:,:) :: gleafmasrow
    real, pointer, dimension(:,:,:) :: bleafmasrow
    real, pointer, dimension(:,:,:) :: stemmassrow
    real, pointer, dimension(:,:,:) :: rootmassrow
    real, pointer, dimension(:,:) :: twarmm
    real, pointer, dimension(:,:) :: tcoldm
    real, pointer, dimension(:,:) :: gdd5
    real, pointer, dimension(:,:) :: aridity           !< aridity index, ratio of potential evaporation to precipitation
    real, pointer, dimension(:,:) :: srplsmon          !< number of months in a year with surplus water i.e.precipitation more than potential evaporation
    real, pointer, dimension(:,:) :: defctmon          !< number of months in a year with water deficit i.e.precipitation less than potential evaporation
    real, pointer, dimension(:,:) :: anndefct          !< annual water deficit (mm)
    real, pointer, dimension(:,:) :: annsrpls          !< annual water surplus (mm)
    real, pointer, dimension(:,:) :: annpcp            !< annual precipitation (mm)
    real, pointer, dimension(:,:) :: dry_season_length !< length of dry season (months)
    real, pointer, dimension(:,:,:) :: litrmassrow
    real, pointer, dimension(:,:,:) :: soilcmasrow
    integer, pointer, dimension(:,:,:) :: lfstatusrow
    integer, pointer, dimension(:,:,:) :: pandaysrow
    real, pointer, dimension(:,:) :: Cmossmas          !< C in moss biomass, \f$kg C/m^2\f$
    real, pointer, dimension(:,:) :: litrmsmoss        !< moss litter mass, \f$kg C/m^2\f$
    real, pointer, dimension(:,:) :: dmoss             !< depth of living moss (m)


    ! local variables for this subroutine
    integer                           :: checksum, k !, PFTCover_size
    real, allocatable, dimension(:)   :: PFTCover, Soil, Canopy_real, Snow, CPools, Peatlands, CompetClimate !< Checksum values to track
    integer, allocatable, dimension(:) :: Canopy_int
    character(len = 10)                 :: lonchar, latchar, checksumchar !< string representations
    character(len = 500)                :: buffer, filename

    ! point pointers:
    ctem_on           => c_switch%ctem_on
    PFTCompetition    => c_switch%PFTCompetition
    lnduseon          => c_switch%lnduseon
    fcancmxrow        => vrot%fcancmx
    gleafmasrow       => vrot%gleafmas
    bleafmasrow       => vrot%bleafmas
    stemmassrow       => vrot%stemmass
    rootmassrow       => vrot%rootmass
    twarmm            => vrot%twarmm
    tcoldm            => vrot%tcoldm
    gdd5              => vrot%gdd5
    aridity           => vrot%aridity
    srplsmon          => vrot%srplsmon
    defctmon          => vrot%defctmon
    anndefct          => vrot%anndefct
    annsrpls          => vrot%annsrpls
    annpcp            => vrot%annpcp
    dry_season_length => vrot%dry_season_length
    litrmassrow       => vrot%litrmass
    soilcmasrow       => vrot%soilcmas
    lfstatusrow       => vrot%lfstatus
    pandaysrow        => vrot%pandays
    Cmossmas          => vrot%Cmossmas
    litrmsmoss        => vrot%litrmsmoss
    dmoss             => vrot%dmoss
    FCANROT           => class_rot%FCANROT
    FAREROT           => class_rot%FAREROT
    TBARROT           => class_rot%TBARROT
    THLQROT           => class_rot%THLQROT
    THICROT           => class_rot%THICROT
    TCANROT           => class_rot%TCANROT
    TSNOROT           => class_rot%TSNOROT
    TPNDROT           => class_rot%TPNDROT
    ZPNDROT           => class_rot%ZPNDROT
    RCANROT           => class_rot%RCANROT
    SCANROT           => class_rot%SCANROT
    SNOROT            => class_rot%SNOROT
    ALBSROT           => class_rot%ALBSROT
    RHOSROT           => class_rot%RHOSROT
    GROROT            => class_rot%GROROT


    write(lonchar, '(I4)')lonIndex
    write(latchar, '(I4)')latIndex
    ! generate the proper, formatted filename
    !ignoreLint(3) (it messes with the file path)
    write(filename, "(A,A,A,'_',A,A)") TRIM(adjustl(c_switch%output_directory)), '/checksums/', &
                                       TRIM(adjustl(lonchar)), TRIM(adjustl(latchar)), '.csv'
    open(unit = 500, file = TRIM(adjustl(filename)))


    checksum = 0
    allocate(PFTCover(SIZE(FAREROT) + SIZE(FCANROT) + SIZE(GROROT) + SIZE(fcancmxrow))) !< allocate the array
    PFTCover = (/ pack(FAREROT, .true.), pack(FCANROT, .true.), pack(GROROT, .true.), pack(fcancmxrow, .true.) /)
    do k = 1, SIZE(PFTCover)
      checksum = checksum + bitcount(PFTCover(k))
    end do
    write(checksumchar, '(I4)')checksum   !< transfer to a string variable
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)), TRIM(adjustl(latchar)), &
                                          "PFTCover",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(Soil(SIZE(THLQROT) + SIZE(THICROT) + SIZE(TPNDROT) + SIZE(ZPNDROT)))
    Soil = (/ pack(THLQROT, .true.), pack(THICROT, .true.), pack(TPNDROT, .true.), pack(ZPNDROT, .true.) /)
    do k = 1, SIZE(Soil)
      checksum = checksum + bitcount(Soil(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)), TRIM(adjustl(latchar)), &
                                          "Soil",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(Canopy_real(SIZE(TCANROT) + SIZE(RCANROT) + SIZE(SCANROT)))
    Canopy_real = (/ pack(TCANROT, .true.), pack(RCANROT, .true.), pack(SCANROT, .true.) /)
    allocate(Canopy_int(SIZE(lfstatusrow) + SIZE(pandaysrow)))
    Canopy_int = (/ pack(lfstatusrow, .true.), pack(pandaysrow, .true.) /)
    do k = 1, SIZE(Canopy_real)
      checksum = checksum + bitcount(Canopy_real(k))
    end do
    do k = 1, SIZE(Canopy_int)
      checksum = checksum + bitcount_int(Canopy_int(k))
    end do

    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)), TRIM(adjustl(latchar)), &
                                          "Canopy",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(Snow(SIZE(TSNOROT) + SIZE(SNOROT) +  SIZE(ALBSROT) + SIZE(RHOSROT)))
    Snow = (/ pack(TSNOROT, .true.), pack(SNOROT, .true.), pack(ALBSROT, .true.), pack(RHOSROT, .true.) /)
    do k = 1, SIZE(Snow)
      checksum = checksum + bitcount(Snow(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)), TRIM(adjustl(latchar)), &
                                          "Snow",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(CPools(SIZE(gleafmasrow) + SIZE(bleafmasrow) + SIZE(stemmassrow) + SIZE(rootmassrow) + SIZE(litrmassrow) + SIZE(soilcmasrow)))
    CPools = (/ pack(gleafmasrow, .true.), pack(bleafmasrow, .true.), pack(stemmassrow, .true.), pack(rootmassrow, .true.), pack(litrmassrow, .true.), pack(soilcmasrow, .true.) /)
    do k = 1, SIZE(Cpools)
      checksum = checksum + bitcount(Cpools(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)), TRIM(adjustl(latchar)), &
                                          "Cpools",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(Peatlands(SIZE(Cmossmas) + SIZE(litrmsmoss) + SIZE (dmoss)))
    Peatlands = (/ pack(Cmossmas, .true.), pack(litrmsmoss, .true.), pack(dmoss, .true.) /)
    do k = 1, SIZE(Peatlands)
      checksum = checksum + bitcount(Peatlands(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)), TRIM(adjustl(latchar)), &
                                          "Peatlands",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(CompetClimate(SIZE(twarmm) + SIZE(tcoldm) + SIZE(gdd5) + SIZE(aridity) + SIZE(srplsmon) + SIZE(defctmon) + SIZE(anndefct) + SIZE(annsrpls) + SIZE(annpcp) + SIZE(dry_season_length)))
    CompetClimate = (/ pack(twarmm, .true.), pack(tcoldm, .true.), pack(gdd5, .true.), pack(aridity, .true.), pack(srplsmon, .true.), pack(defctmon, .true.), pack(anndefct, .true.), pack(annsrpls, .true.), pack(annpcp, .true.), pack(dry_season_length, .true.) /)
    do k = 1, SIZE(CompetClimate)
      checksum = checksum + bitcount(CompetClimate(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)), TRIM(adjustl(latchar)), &
                                          "CompetClimate",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    close(500)

  end subroutine checksumCalc
  !! @}

  !> \ingroup generalUtils_bitcount
  !! @{
  !! This function generates the bitcount of the specified variable
  !! @author M. Fortier
  !!
  integer function bitcount(scalar)

    implicit none
    real    :: scalar     !< Scalar to be bit-counted
    integer(SELECTED_REAL_KIND(15, 307)) :: scalar_int !< Integer with memory representation of 'scalar'

    integer :: bit
    bitcount = 0
    scalar_int = TRANSFER(scalar, scalar_int)
    do bit = 0, BIT_SIZE(scalar_int) - 1
      if (BTEST(scalar_int, bit) ) bitcount = bitcount + 1
    end do

  end function bitcount
  !! @}


  !> \ingroup generalUtils_bitcount_int
  !! @{
  !! This function generates the bitcount of the specified variable, from an integer
  !! @author M. Fortier
  !!
  integer function bitcount_int(scalar)

    implicit none
    integer    :: scalar     !< Scalar to be bit-counted
    integer    :: bit

    bitcount_int = 0
    do bit = 0, BIT_SIZE(scalar) - 1
      if (BTEST(scalar, bit) ) bitcount_int = bitcount_int + 1
    end do

  end function bitcount_int
  !! @}
  ! -----------------------------------------------------------------------------------------------
  !> \namespace checksum

  !>
  !! @author M. Fortier
  !!
  !! This module computes content-based checksums of several groups of variables after
  !! a run has completed. These checksums are used when making non-logical changes to
  !! the code-base. This is accomplished by comparing the computed checksums, to the
  !! checksums of an identical run (same parameters, met files, etc.) before the changes
  !! were made. This is a modified form of regression testing specific to non-logical
  !! software changes.
  !!
  !! Checksums are not an infallible means of ensuring no logical changes, as two numbers
  !! may have different binary representations with the same number of flipped bits. However,
  !! with the number of variables we are checking, it is extremely unlikely to render a false-
  !! positive. If we take a data item of bit-length \f$n\f$, having \f$b\f$ flipped bits in
  !! its representation, the number of same-length data items with the same checksum can
  !! be expressed through the binary coefficient
  !!
  !! \f[ {n\choose b} \f]
  !!
  !! If we assume the worst case where \f$b=\frac{n}{2}\f$, then we are left with
  !!
  !! \f[ {n\choose n/2} \f]
  !!
  !! Dividing by the total number of possible values for an n-digit binary number,
  !! we find the probability of a false positive checksum to be, in the worst case:
  !!
  !! \f[ \large \frac{{n\choose n/2}}{2^n} \f]
  !!
  !! Due to the cascading effects of any logical changes to the model, all 7 groups
  !! of variables (6 if peatlands is diabled) will be affected. The smallest of these
  !! groups consists of 3 arrays of numbers for which checksums are computed. With the
  !! conservative assumption that each array contains only a single 32-bit number, the
  !! probability of a false positive (in the absolute worst case) becomes:
  !!
  !! \f[ \large \left({\frac{96\choose 48}{2^{96}}}\right)^{6} ~= 2.87105123*10^{-7}       \f]
  !!
  !> \file
end module checksum
