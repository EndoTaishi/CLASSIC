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
  subroutine tracerDynamics()
    
    use ctem_statevars, only : c_switch

    implicit none 
    
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
      call simpleTracer()
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
  subroutine simpleTracer()
    
    use classic_params, only : icc, deltat
    use ctem_statevars, only : tracer,c_switch,vgat
    
    implicit none 
    

    real, pointer :: tracerGLeafMass(:,:)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerBLeafMass(:,:)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerStemMass(:,:)       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerRootMass(:,:)       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer :: tracerMossCMass(:)      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, pointer :: tracerMossLitrMass(:)   !< Tracer mass in moss litter, \f$kg C/m^2\f$

    real, pointer :: tracerco2conc(:)
    
    real, pointer :: rmlveg(:,:)      !< Leaf maintenance resp. rate for each pft
    real, pointer :: gppveg(:,:)      !< Gross primary productity for each pft
    real, pointer :: rmsveg(:,:)      !< Stem maintenance resp. rate for each pft
    real, pointer :: rmrveg(:,:)      !< Root maintenance resp. rate for each pft
    real, pointer :: rgveg(:,:)       !< Growth resp. rate for each pft
    real, pointer :: litresveg(:,:,:)  !<fluxes for each pft: litter respiration for each pft + bare fraction
    real, pointer :: soilcresveg(:,:,:)  !<soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts
    real, pointer :: humiftrsveg(:,:,:) !<
    real, pointer :: reprocost(:,:)   !< Cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) 
    real, pointer :: afrleaf(:,:)     !<allocation fraction for leaves
    real, pointer :: afrstem(:,:)     !<allocation fraction for stems
    real, pointer :: afrroot(:,:)     !<allocation fraction for roots
    real, pointer :: tltrleaf(:,:)    !<total leaf litter fall rate (u-mol co2/m2.sec)
    real, pointer :: tltrstem(:,:)    !<total stem litter fall rate (u-mol co2/m2.sec)
    real, pointer :: tltrroot(:,:)    !<total root litter fall rate (u-mol co2/m2.sec)

    ! Local
    
    ! Point pointers 

    tracerco2conc     => tracer%tracerCO2rot
    tracerGLeafMass   => tracer%gLeafMassrow
    tracerBLeafMass   => tracer%bLeafMassrow
    tracerStemMass    => tracer%stemMassrow
    tracerRootMass    => tracer%rootMassrow
    tracerLitrMass    => tracer%litrMassrow
    tracerSoilCMass   => tracer%soilCMassrow
    tracerMossCMass   => tracer%mossCMassrow
    tracerMossLitrMass => tracer%mossLitrMassrow
    
    rmlveg            => vgat%rmlvegacc
    rmrveg            => vgat%rmrveg
    rmsveg            => vgat%rmsveg
    gppveg            => vgat%gppveg
    rgveg             => vgat%rgveg
    litresveg         => vgat%litresveg
    soilcresveg       => vgat%soilcresveg
    humiftrsveg       => vgat%humiftrsveg
    reprocost         => vgat%reprocost
    afrleaf           => vgat%afrleaf
    afrstem           => vgat%afrstem
    afrroot           => vgat%afrroot
    tltrleaf          => vgat%tltrleaf
    tltrstem          => vgat%tltrstem
    tltrroot          => vgat%tltrroot
    
    
    ! ---------
    
    
    ! u mol co2/m2.sec -> \f$(kg C/m^2)\f$ respired over the model time step
    !ltrestep(i,j,k) = ltresveg(i,j,k) * deltat / 963.62
    !screstep(i,j,k) = scresveg(i,j,k) * deltat / 963.62

    
    ! First update the green leaves
    ! tracerGLeafMass + GPP - allocation to roots - allocation to stem - respiration - litterfall - conversion to brown leaves - combusted - growth resp. - repro cost.
    
    ! Update brown leaves 
    ! tracerBLeafMass - respiration - litterfall - combusted
    
    ! Update stem mass 
    ! tracerStemMass + allocation - respiration - litter fall  - combusted
    
    ! Update root mass
    ! tracerRootMass + allocation - respiration - litter fall - combusted
    
    ! Update litter mass 
    ! tracerLitrMass + litterfall (leaves, stem, roots) - respiration - humification
    
    ! Update soil C mass 
    ! tracerSoilCMass + humication - respiration 

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
