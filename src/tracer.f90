!>\file
!> Module containing all relevant subroutines for the model 
!! tracer.
!!@author Joe Melton
!!
module tracerModule

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
!!
!! The tracer trackes the C movement through the green leaf,
!! brown leaf, root, stem, litter and soil C. Carbon that is 
!! incorporated into the plants are given a tracer value of 
!! tracerValue, which corresponds to that read in from the 
!! tracerCO2 file for the year simulated. As the simulation 
!! runs and C is transferred from the living pools to the detrital 
!! pools, the tracer also is transfered. 
!> @author Joe Melton
  subroutine simpleTracer(il1,il2)
    
    use classic_params, only : icc, deltat, iccp2, zero, grass, ignd
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
    
    real, pointer :: litresveg(:,:,:)  !<fluxes for each pft: litter respiration for each pft + bare fraction
    real, pointer :: soilcresveg(:,:,:)  !<soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts
    real, pointer :: humiftrsveg(:,:,:) !<
    real, pointer :: reprocost(:,:)   !< Cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) 
    real, pointer :: tltrleaf(:,:)    !<total leaf litter fall rate (u-mol co2/m2.sec)
    real, pointer :: tltrstem(:,:)    !<total stem litter fall rate (u-mol co2/m2.sec)
    real, pointer :: tltrroot(:,:)    !<total root litter fall rate (u-mol co2/m2.sec)
    real, pointer :: blfltrdt(:,:)    !<brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, pointer :: glfltrdt(:,:)    !<green leaf litter generated due to disturbance \f$(kg c/m^2)\f$
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
    real, pointer :: mortLeafGtoB(:,:)  !< Green leaf mass converted to brown due to mortality \f$(kg C/m^2)\f$
    real, pointer :: phenLeafGtoB(:,:)  !< Green leaf mass converted to brown due to phenology \f$(kg C/m^2)\f$
    real, pointer :: fcancmx(:,:)    !< Maximum fractional coverage of CTEM PFTs, but this can be
                                     !! modified by land-use change, and competition between PFTs
    real, pointer :: leaflitr(:,:)   !< Leaf litter fall rate (\f$\mu mol CO2 m^{-2} s^{-1}\f$). 
                                     !! this leaf litter does not include litter generated 
                                     !! due to mortality/fire
    real, pointer :: rmatctem(:,:,:)   !< Fraction of roots for each of CTEM's PFTs in each soil layer
    real, pointer :: turbLitter(:,:,:) !< Litter gains/losses due to turbation [ \f$kg C/m^2\f$ ], negative is a gain.
    real, pointer :: turbSoilC(:,:,:)  !< Soil C gains/losses due to turbation [ \f$kg C/m^2\f$ ], negative is a gain.
    real, pointer :: gleafmas(:,:)     !< Green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, pointer :: bleafmas(:,:)     !< Brown leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, pointer :: stemmass(:,:)     !< Stem mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, pointer :: rootmass(:,:)     !< Root mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, pointer :: litrmass(:,:,:)   !< Litter mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    real, pointer :: soilcmas(:,:,:)   !< Soil carbon mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    real, pointer, dimension(:,:) :: gLeafLandCompChg    !< Tracker variable for C movement due to competition and LUC in the green leaf pool  [ \f$kg C/m^2\f$ ], negative is a gain.
    real, pointer, dimension(:,:) :: bLeafLandCompChg    !< Tracker variable for C movement due to competition and LUC in the brown leaf pool  [ \f$kg C/m^2\f$ ], negative is a gain.
    real, pointer, dimension(:,:) :: stemLandCompChg   !< Tracker variable for C movement due to competition and LUC in the stem pool  [ \f$kg C/m^2\f$ ], negative is a gain.
    real, pointer, dimension(:,:) :: rootLandCompChg   !< Tracker variable for C movement due to competition and LUC in the root pool  [ \f$kg C/m^2\f$ ], negative is a gain.
    real, pointer, dimension(:,:,:) :: litterLandCompChg !< Tracker variable for C movement due to competition and LUC in the litter pool  [ \f$kg C/m^2\f$ ], negative is a gain.
    real, pointer, dimension(:,:,:) :: soilCLandCompChg !< Tracker variable for C movement due to competition and LUC in the soil C pool  [ \f$kg C/m^2\f$ ], negative is a gain.

    
    ! Local
    integer :: i,j,k
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
    litresveg         => vgat%litresveg
    soilcresveg       => vgat%soilcresveg
    humiftrsveg       => vgat%humiftrsveg
    reprocost         => vgat%reprocost
    tltrleaf          => vgat%tltrleaf
    tltrstem          => vgat%tltrstem
    tltrroot          => vgat%tltrroot
    blfltrdt          => vgat%blfltrdt
    glfltrdt          => vgat%glfltrdt
    glcaemls          => vgat%glcaemls
    blcaemls          => vgat%blcaemls
    rtcaemls          => vgat%rtcaemls
    stcaemls          => vgat%stcaemls
    ltrcemls          => vgat%ltrcemls
    ntchlveg          => vgat%ntchlveg
    ntchsveg          => vgat%ntchsveg
    ntchrveg          => vgat%ntchrveg
    mortLeafGtoB      => vgat%mortLeafGtoB
    phenLeafGtoB      => vgat%phenLeafGtoB
    fcancmx           => vgat%fcancmx
    leaflitr          => vgat%leaflitr
    rmatctem          => vgat%rmatctem
    turbLitter        => vgat%turbLitter
    turbSoilC         => vgat%turbSoilC
    gleafmas          => vgat%gleafmas
    bleafmas          => vgat%bleafmas
    stemmass          => vgat%stemmass
    rootmass          => vgat%rootmass
    litrmass          => vgat%litrmass
    soilcmas          => vgat%soilcmas
    gLeafLandCompChg  => vgat%gLeafLandCompChg
    bLeafLandCompChg  => vgat%bLeafLandCompChg
    stemLandCompChg  => vgat%stemLandCompChg
    rootLandCompChg  => vgat%rootLandCompChg
    litterLandCompChg => vgat%litterLandCompChg
    soilCLandCompChg  => vgat%soilCLandCompChg
    ! ---------
    
    convertUnits = deltat / 963.62
    
    !> Since this is a simple tracer, we don't need to do any conversion of the 
    !! carbon that is uptaked in this timestep. We assume that the tracerValue
    !! should just be applied as is.
    do i = il1, il2 
      do j = 1, iccp2
        if (j <= icc) then !these are just icc sized arrays.
          if (fcancmx(i,j) > zero) then
                        
            ! *** Update the green leaves *** 
            
            ! FLAG at present I don't account for LUC!
            
            if (ntchlveg(i,j) > 0.) then ! NPP was positive to leaves
              gains = ntchlveg(i,j) * convertUnits * 1. !tracerValue(i)  
            else ! loss of C from leaves due to negative NPP.
              gains = ntchlveg(i,j) * convertUnits  
            end if
            ! tltrleaf (phenology litter, mortality, disturbance) includes brown leaf
            ! generation from disturbance so remove that)
            if (.not. grass(j) ) then
              losses = (tltrleaf(i,j) - (blfltrdt(i,j)/convertUnits)) * convertUnits & 
                       + glcaemls(i,j) * convertUnits !combusted by fire.
            else 
              losses = glfltrdt(i,j) & !green leaf litter generated by fire.
                       + mortLeafGtoB(i,j) & !mortality transfer to brown leaves, only >0 for grasses
                       + phenLeafGtoB(i,j) & !phenology transfer to brown leaves, only >0 for grasses
                       + glcaemls(i,j) * convertUnits !combusted by fire.
            end if 
                     
            !> When the tracer is calculated we apply the gain and losses to the exisiting
            !! pool. We also include the gains/losses due to changing PFT areas after land 
            !! use change or competition. 
            
            !! FLAG OK to assign no special tracer value to LUC/Comp???? *********
            tracerGLeafMass(i,j) = tracerGLeafMass(i,j) + gains - losses - gLeafLandCompChg(i,j)
            if (tracerGLeafMass(i,j) < zero) tracerGLeafMass(i,j) = 0.
            
            ! ***  Update brown leaves (grass only) *** 
            if (grass(j) ) then
              gains = phenLeafGtoB(i,j) & ! phenology transfer to brown leaves, only >0 for grasses
                      + mortLeafGtoB(i,j) ! mortality transfer to brown leaves, only >0 for grasses
              
              losses = leaflitr(i,j)  * convertUnits & ! phenology litter generation.
                       + blfltrdt(i,j) & !brown leaf litter generated by fire.
                       + blcaemls(i,j) * convertUnits !combusted by fire.
              
              tracerBLeafMass(i,j) = tracerBLeafMass(i,j) + gains - losses - bLeafLandCompChg(i,j)
              if (tracerBLeafMass(i,j) < zero) tracerBLeafMass(i,j) = 0.
            end if 
            
            ! ***  Update stem mass (all except grass) *** 
            if (.not. grass(j)) then
            
              if (ntchsveg(i,j) > 0.) then ! NPP was positive to leaves
                gains = ntchsveg(i,j) * convertUnits * 1. !tracerValue(i)  
              else ! loss of C from leaves due to negative NPP.
                gains = ntchsveg(i,j) * convertUnits  
              end if

              losses = tltrstem(i,j) * convertUnits & !turnover, mortality, disturbance => litter.
                      + stcaemls(i,j) * convertUnits ! combusted by fire.
              
              tracerStemMass(i,j) = tracerStemMass(i,j) + gains - losses - stemLandCompChg(i,j)
              if (tracerStemMass(i,j) < zero) tracerStemMass(i,j) = 0.
              
            end if 
            
            ! ***  Update root mass *** 
            if (ntchrveg(i,j) > 0.) then ! NPP was positive to leaves
              gains = ntchrveg(i,j) * convertUnits * 1. !tracerValue(i)  
            else ! loss of C from leaves due to negative NPP.
              gains = ntchrveg(i,j) * convertUnits  
            end if
            
            losses = tltrroot(i,j) * convertUnits &
                    + rtcaemls(i,j) * convertUnits
                    
            tracerRootMass(i,j) = tracerRootMass(i,j) + gains - losses - rootLandCompChg(i,j)
            if (tracerRootMass(i,j) < zero) tracerRootMass(i,j) = 0.

          else 
            ! Set the tracer pools to 0 if the model pool is zero
            ! or the PFT is not in the tile.          
            tracerLitrMass(i,j,:) = 0.
            tracerSoilCMass(i,j,:) = 0.
            tracerRootMass(i,j) = 0.
            tracerStemMass(i,j) = 0.
            tracerBLeafMass(i,j) = 0.
            tracerGLeafMass(i,j) = 0.
          end if         
          if (rootmass(i,j) < zero) tracerRootMass(i,j) = 0.
          if (stemmass(i,j) < zero) tracerStemMass(i,j) = 0.
          if (bleafmas(i,j) < zero) tracerBLeafMass(i,j) = 0.
          if (gleafmas(i,j) < zero) tracerGLeafMass(i,j) = 0.
        end if
      
        ! ***  Update litter mass *** 
        do k = 1, ignd
          if (j <= icc) then 
            if (k == 1) then 
              ! surface gains 
              gains = (tltrleaf(i,j) & !leaf litter 
                      + tltrstem(i,j) & ! stem 
                      + tltrroot(i,j) * rmatctem(i,j,k) & ! root 
                      + reprocost(i,j))  * convertUnits  ! reproductive tissues
              losses = (litresveg(i,j,k)  & !litter respiration
                      + humiftrsveg(i,j,k) &  ! humification
                      + ltrcemls(i,j)) * convertUnits !combusted
            else 
              ! deeper layers 
              gains = tltrroot(i,j) * rmatctem(i,j,k)  * convertUnits !from roots only.
              losses = (litresveg(i,j,k)  & !litter respiration
                      + humiftrsveg(i,j,k)) * convertUnits !humification
            end if                
          else ! bareground or LUC product pools
            gains = 0. ! The gains come from LUC, so are in the litterLandCompChg term.
            losses = (litresveg(i,j,k)  & !litter respiration
                    + humiftrsveg(i,j,k)) * convertUnits  ! humification
          end if 
          
          ! Litter mass is mixed by soil turbation (cryo, bio) so we also consider the turbation movements
          ! but LUC product pools are not considered to be turbated.
          if (j < iccp2) then
            tracerLitrMass(i,j,k) = tracerLitrMass(i,j,k) + gains - losses &
                                                          - turbLitter(i,j,k) &
                                                          - litterLandCompChg(i,j,k)
          else 
            tracerLitrMass(i,j,k) = tracerLitrMass(i,j,k) + gains - losses  - litterLandCompChg(i,j,k)
          end if 
          if (tracerLitrMass(i,j,k) < zero) tracerLitrMass(i,j,k) = 0.
        end do 
                
        ! ***  Update soil C mass ***          
        do k = 1, ignd
          gains = humiftrsveg(i,j,k) * convertUnits ! humification 
          losses = soilcresveg(i,j,k) * convertUnits ! respiration 
          if (j < iccp2) then
            tracerSoilCMass(i,j,k) = tracerSoilCMass(i,j,k) + gains - losses &
                                                            - turbSoilC(i,j,k)&
                                                            - soilCLandCompChg(i,j,k)         
          else ! no turbation of LUC product pools.
            tracerSoilCMass(i,j,k) = tracerSoilCMass(i,j,k) + gains - losses - soilCLandCompChg(i,j,k)
          end if
          if (tracerSoilCMass(i,j,k) < zero) tracerSoilCMass(i,j,k) = 0.
        end do
        
        if (sum(litrmass(i,j,:)) < zero) tracerLitrMass(i,j,:) = 0.
        if (sum(soilcmas(i,j,:)) < zero) tracerSoilCMass(i,j,:) = 0.       
           
      end do ! j 
    end do ! i 
    
    ! k=1
    ! write(*,'(a7,i4,2f13.10)')'tracer1',k,tracerLitrMass(1,6,k),tracerLitrMass(1,7,k)
    ! write(*,'(a7,i4,2f13.10)')'litter1',k,litrmass(1,6,k),litrmass(1,7,k)      
    ! k=2  
    ! write(*,'(a7,i4,2f13.10)')'tracer2',k,tracerLitrMass(1,6,k),tracerLitrMass(1,7,k)
    ! write(*,'(a7,i4,2f13.10)')'litter2',k,litrmass(1,6,k),litrmass(1,7,k)    
    ! k=3    
    ! write(*,'(a7,i4,2f13.10)')'tracer3',k,tracerLitrMass(1,6,k),tracerLitrMass(1,7,k)
    ! write(*,'(a7,i4,2f13.10)')'litter3',k,litrmass(1,6,k),litrmass(1,7,k)  
    ! k=4     
    ! write(*,'(a7,i4,2f13.10)')'tracer4',k,tracerLitrMass(1,6,k),tracerLitrMass(1,7,k)
    ! write(*,'(a7,i4,2f13.10)')'litter4',k,litrmass(1,6,k),litrmass(1,7,k)        

    !read(*,*)
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
end module tracerModule
