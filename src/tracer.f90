!>\file
!> Module containing all relevant subroutines for the model 
!! tracer.
!!@author Joe Melton
!!
module tracer

implicit none

public :: tracerDynamics
public :: simpleTracer
public :: tracer13C
public :: tracer14C
!public :: tracerBalance

contains

!>\ingroup tracer_tracerDynamics
!!@{Determine which tracer subroutine is to be called based on 
!! value of useTracer. Either simpleTracer, 13C, or 14C.
!> @author Joe Melton
  subroutine tracerDynamics(il1,il2)
    
    use ctem_statevars, only : c_switch

    implicit none 
    
    integer, intent(in) :: il1 !<other variables: il1=1
    integer, intent(in) :: il2 !<other variables: il2=ilg

    integer, pointer :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used. 
                                  !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
                                  !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.                         
                                  !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme. 
                                  !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.  
    useTracer         => c_switch%useTracer
    
    ! Determine which tracer subroutine we are interested in.
    if (useTracer == 0) then 
      print*,'Error, entered tracerDynamics with a useTracer value of 0'
      call xit('tracer',0)
    else if (useTracer == 1) then 
      call simpleTracer(il1,il2)
    else if (useTracer == 2) then 
      call tracer13C()
    else if (useTracer == 3) then 
      call tracer14C()
    end if
    
  end subroutine tracerDynamics
!!@}
! -------------------------------------------------------
!>\ingroup tracer_simpleTracer
!!@{ Simple tracer which tracks C flow through the system.
!! No fractionation effects. The tracer's value depends on how 
!! the model is initialized and the input file used.
!> @author Joe Melton
  subroutine simpleTracer(il1,il2)
    
    use classic_params, only : icc, deltat, iccp2
    use ctem_statevars, only : tracer,c_switch,vgat
    
    implicit none 
    
    integer, intent(in) :: il1 !<other variables: il1=1
    integer, intent(in) :: il2 !<other variables: il2=ilg

    real, pointer :: tracerGLeafMass(:,:)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerBLeafMass(:,:)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerStemMass(:,:)       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerRootMass(:,:)       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer :: tracerMossCMass(:)      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, pointer :: tracerMossLitrMass(:)   !< Tracer mass in moss litter, \f$kg C/m^2\f$

    real, pointer :: tracerValue(:)
    
    ! real, pointer :: rmlveg(:,:)      !< Leaf maintenance resp. rate for each pft
    ! real, pointer :: gppveg(:,:)      !< Gross primary productity for each pft
    ! real, pointer :: rmsveg(:,:)      !< Stem maintenance resp. rate for each pft
    ! real, pointer :: rmrveg(:,:)      !< Root maintenance resp. rate for each pft
    ! real, pointer :: rgveg(:,:)       !< Growth resp. rate for each pft
    real, pointer :: litresveg(:,:,:)  !<fluxes for each pft: litter respiration for each pft + bare fraction
    real, pointer :: soilcresveg(:,:,:)  !<soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts
    real, pointer :: humiftrsveg(:,:,:) !<
    real, pointer :: reprocost(:,:)   !< Cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) 
    ! real, pointer :: afrleaf(:,:)     !<allocation fraction for leaves
    ! real, pointer :: afrstem(:,:)     !<allocation fraction for stems
    ! real, pointer :: afrroot(:,:)     !<allocation fraction for roots
    real, pointer :: tltrleaf(:,:)    !<total leaf litter fall rate (u-mol co2/m2.sec)
    real, pointer :: tltrstem(:,:)    !<total stem litter fall rate (u-mol co2/m2.sec)
    real, pointer :: tltrroot(:,:)    !<total root litter fall rate (u-mol co2/m2.sec)
    real, pointer :: blfltrdt(:,:)    !<brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, pointer :: glcaemls(:,:)  !<green leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, pointer :: blcaemls(:,:)  !<brown leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, pointer :: rtcaemls(:,:)  !<root carbon emission disturbance losses, \f$kg c/m^2\f$
    real, pointer :: stcaemls(:,:)  !<stem carbon emission disturbance losses, \f$kg c/m^2\f$
    real, pointer :: ltrcemls(:,:)  !<litter carbon emission disturbance losses, \f$kg c/m^2\f$
    real, pointer :: ntchlveg(:,:)  !<fluxes for each pft: Net change in leaf biomass, u-mol CO2/m2.sec
    real, pointer :: ntchsveg(:,:)  !<fluxes for each pft: Net change in stem biomass, u-mol CO2/m2.sec
    real, pointer :: ntchrveg(:,:)  !<fluxes for each pft: Net change in root biomass, 
                                    !! the net change is the difference between allocation and
                                    !! autotrophic respiratory fluxes, u-mol CO2/m2.sec

    ! Local
    integer :: i,j
    real :: gains, losses
    real :: convertUnits      !< This converts the units from u-mol CO2/m2.sec to kg C/m^2
    
    ! Point pointers 
    tracerValue       => tracer%tracerCO2rot
    tracerGLeafMass   => tracer%gLeafMassgat
    tracerBLeafMass   => tracer%bLeafMassgat
    tracerStemMass    => tracer%stemMassgat
    tracerRootMass    => tracer%rootMassgat
    tracerLitrMass    => tracer%litrMassgat
    tracerSoilCMass   => tracer%soilCMassgat
    tracerMossCMass   => tracer%mossCMassgat
    tracerMossLitrMass => tracer%mossLitrMassgat
    ! rmlveg            => vgat%rmlvegacc
    ! rmrveg            => vgat%rmrveg
    ! rmsveg            => vgat%rmsveg
    ! gppveg            => vgat%gppveg
    ! rgveg             => vgat%rgveg
    litresveg         => vgat%litresveg
    soilcresveg       => vgat%soilcresveg
    humiftrsveg       => vgat%humiftrsveg
    reprocost         => vgat%reprocost
    ! afrleaf           => vgat%afrleaf
    ! afrstem           => vgat%afrstem
    ! afrroot           => vgat%afrroot
    tltrleaf          => vgat%tltrleaf
    tltrstem          => vgat%tltrstem
    tltrroot          => vgat%tltrroot
    blfltrdt          => vgat%blfltrdt
    glcaemls          => vgat%glcaemls
    blcaemls          => vgat%blcaemls
    rtcaemls          => vgat%rtcaemls
    stcaemls          => vgat%stcaemls
    ltrcemls          => vgat%ltrcemls
    ntchlveg          => vgat%ntchlveg
    ntchsveg          => vgat%ntchsveg
    ntchrveg          => vgat%ntchrveg

    ! ---------
    
    convertUnits = deltat / 963.62
    
    !> Since this is a simple tracer, we don't need to do any conversion of the 
    !! carbon that is uptaked in this timestep. We assume that the tracerValue
    !! should just be applied as is.
    do i = il1, il2 
      do j = 1, iccp2
        if (j <= icc) then !these are just icc sized arrays.

          ! First update the green leaves
          if (ntchlveg(i,j) > 0.) then ! NPP was positive to leaves
            gains = ntchlveg(i,j) * tracerValue(i)  * convertUnits 
          else ! loss of C from leaves due to negative NPP.
            gains = ntchlveg(i,j) * convertUnits  
          end if
          ! Losses include: leaf loss to mortality, phenology, and fire (both conversion to litter and combustion)
          ! tltrleaf includes brown leaf generation from disturbance so remove that)
          losses = (tltrleaf(i,j) - (blfltrdt(i,j)/convertUnits) & 
                + glcaemls(i,j)) * convertUnits 
                
          tracerGLeafMass(i,j) = tracerGLeafMass(i,j) + gains - losses
          
          ! Update brown leaves 
          ! tracerBLeafMass - respiration - litterfall - combusted
          
          ! Update stem mass 
          ! tracerStemMass + allocation - respiration - litter fall  - combusted
          
          ! Update root mass
          ! tracerRootMass + allocation - respiration - litter fall - combusted
      
        end if
      
      ! Update litter mass 
      ! tracerLitrMass + litterfall (leaves, stem, roots) - respiration - humification
      
      ! Update soil C mass 
      ! tracerSoilCMass + humication - respiration 

      end do 
    end do

  end subroutine simpleTracer
!!@}
! -------------------------------------------------------
!>\ingroup tracer_tracer13C
!!@{ 13C tracer -- NOT IMPLEMENTED YET.
  subroutine tracer13C()
    
    use ctem_statevars, only : c_switch

    implicit none 
    
    print*,'tracer13C: Not implemented yet. Usetracer cannot == 2!'
    call xit('tracer',-1)
  end subroutine tracer13C
!!@}
! -------------------------------------------------------
!>\ingroup tracer_tracer13C
!!@{ 13C tracer -- NOT IMPLEMENTED YET.
  subroutine tracer14C()
    
    use ctem_statevars, only : c_switch

    implicit none 
    
    print*,'tracer14C: Not implemented yet. Usetracer cannot == 3!'
    call xit('tracer',-2)
  end subroutine tracer14C
!!@}
! -------------------------------------------------------
!>\namespace tracer
!!
!!
!>\file
end module tracer
