!>\file
!>Applies allometric relationships to converts biomass to 
!! structural attributes
!!@author Vivek Arora, Joe Melton
!!
module applyAllometry

implicit none

public :: allometry
public :: unfrozenRoots

contains

!>\ingroup allometry_bio2str
!!@{
!> Applies allometric relationships to converts biomass to 
!! structural attributes. This subroutine converts leaf, stem, and root biomass 
!! into LAI, vegetation height, and fraction of roots
!! in each soil layer. Storage lai is also calculated.
!!                    
!! Note that while CTEM keeps track of icc pfts, CLASS
!! keeps track of ican pfts, so all these vegetation 
!! structural attributes are calculated for icc pfts
!! separately but then lumped into ican pfts for use in
!! energy & water balance calculations by class
!! also, this subroutine does not estimate zolnc(i,icp1)
!! the log of roughness length over the bare fraction
!! of the grid cell. only roughness lengths over the
!! vegetated fraction are updated

!> @author V. Arora, J. Melton, Y. Peng, S. Sun
  subroutine    allometry( gleafmas, bleafmas, stemmass, rootmass,  & ! In
                             il1,      il2,      ilg,    zbotw,& ! In
                          delzw, soildpth, fcancmx,& ! In
                      ipeatland,maxAnnualActLyr,& ! In
                          ailcg,    ailcb,     ailc,    zolnc,& !Out
                          rmatc, rmatctem,     slai,  bmasveg,& !Out
                       cmasvegc,  veghght, rootdpth,   alvisc,& !Out
                         alnirc,     paic,    slaic) !Out
  !
  !     ----------------------------------------------------------------
  !     16  Aug 2017  - Add BLD COLD-Deciduous shrub 
  !     J. Melton/S. Sun     
  !
  !     2   Jul 2013  - Integreated classic_params module
  !     J. Melton       
  !
  !     22  Nov 2012  - calling this version 1.1 since a fair bit of ctem
  !     V. Arora        subroutines were changed for compatibility with class
  !                     version 3.6 including the capability to run ctem in
  !                     mosaic/tile version along with class.
  !
  !     24  Sep 2012  - add in checks to prevent calculation of non-present
  !     J. Melton       pfts. cleaned up initialization
  !
  !     18  May  2012 - fix bug for soils with the depth close to the
  !     J. Melton       boundary between layers, was resulting in roots 
  !                     being placed incorrectly.
  !
  !     28  Nov. 2011 - make changes for coupling with class 3.5
  !     Y. Peng
  !
  !     31  Aug. 2009 - make changes for coupling with class 3.4. include
  !     V. Arora        new output variable called plant area index
  !                     (paic) and storage leaf area index (slaic) 
  !
  !     14  Mar. 2003 - this subroutine converts leaf, stem, and root biomass 
  !     V. Arora        into lai, vegetation height, and fraction of roots
  !                     in each soil layer. storage lai is also calculated.
  
  use classic_params,        only : ignd, icc, ican, abszero,&
                                l2max,kk, eta, kappa, kn, lfespany, &
                                fracbofg, specsla, abar, avertmas,&
                                alpha, prcnslai, minslai, mxrtdpth,&
                                albvis, albnir,classpfts, nol2pfts, &
                                reindexPFTs        
  use ctemUtilities, only : genSortIndex!, unfrozenRoots
  
  implicit none

  integer ilg !<
  integer il1 !<input: il1=1
  integer il2 !<input: il2=ilg
  integer i, j, k, m, n
  integer sort(icc)
  integer kend

  logical deeproots
    
  real gleafmas(ilg,icc)     !<input: green or live leaf mass in kg c/m2, for the 9 pfts
  real bleafmas(ilg,icc)     !<input: brown or dead leaf mass in kg c/m2, for the 9 pfts
  real stemmass(ilg,icc)     !<input: stem biomass in kg c/m2, for the 9 pfts
  real rootmass(ilg,icc)     !<input: root biomass in kg c/m2, for the 9 pfts
  real ailcg(ilg,icc)        !<output: green lai for ctem's 9 pfts
  real ailcb(ilg,icc)        !<output: brown lai for ctem's 9 pfts. for now we assume only grasses can have brown leaves.
  real ailc(ilg,ican)        !<output: lai to be used by class
  real zolnc(ilg,ican)       !<output: log of roughness length to be used by class
  real paic(ilg, ican)       !<output: plant area index for class' 4 pfts. this is the sum of leaf area index and stem area index.
  real rmatc(ilg,ican,ignd)  !<output: fraction of live roots in each soil layer for each of the class' 4 pfts
  real fcancmx(ilg,icc)      !<input: max. fractional coverages of ctem's 9 pfts. this is different from fcanc and fcancs
                             !<(which may vary with snow depth). fcancmx doesn't change, unless of course its changed by
                             !<land use change or dynamic vegetation.
  real delzw(ilg,ignd)       !<input: thicknesses of the 3 soil layers
  real zbotw(ilg,ignd)       !<input: bottom of soil layers
  real rmatctem(ilg,icc,ignd)!<output: fraction of live roots in each soil layer for each of ctem's 9 pfts
  real slai(ilg,icc)         !<output: storage or imaginary lai for phenology purposes
  real bmasveg(ilg,icc)      !<output: total (gleaf + stem + root) biomass for each ctem pft, kg c/m2
  real cmasvegc(ilg,ican)    !<output: total canopy mass for each of the 4 class pfts. recall that class requires canopy
                             !<mass as an input, and this is now provided by ctem. kg/m2.
  real sai(ilg,icc)          !< 
  real saic(ilg,ican)        !< 
  real sfcancmx(ilg,ican)    !<sum of fcancmxs
  real alvisc(ilg,ican)      !<output: albedo for 4 class pfts simulated by ctem, visible 
  real alnirc(ilg,ican)      !<output: albedo for 4 class pfts simulated by ctem, near ir
  real pai(ilg,icc)          !< 
  real slaic(ilg,ican)       !<output: storage lai. this will be used as min. lai that class sees
                             !<so that it doesn't blow up in its stomatal conductance calculations.
  real sla(icc)              !< 
  real veghght(ilg,icc)      !<output: vegetation height (meters)
  real fcoeff                !< 
  real lnrghlth(ilg,icc)     !< 
  real averough(ilg,ican)    !< 
  real b(icc)                !< 
  real rootdpth(ilg,icc)     !<output: 99% soil rooting depth (meters) both veghght & rootdpth can be used as diagnostics
                             !<to see how vegetation grows above and below ground, respectively.
  real usealpha(ilg,icc)     !< 
  real a(ilg,icc)            !< 
  real useb(ilg,icc)         !< 
  real zroot                 !< 
  real soildpth(ilg)         !<input: soil depth (m)
  real, dimension(ilg), intent(in) :: maxAnnualActLyr   !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).
  real etmp(ilg,icc,ignd)    !<
  real totala(ilg,icc)       !<
  real rmat_sum              !< 

  integer ipeatland(ilg)     !< Peatland flag, non-peatlands are 0.
      
  !! CLASS' original root parameterization has deeper roots than CTEM's
  !! default values based on literature. in the coupled model this leads
  !! to lower evapotranspiration (et) values. an option is provided here to
  !! deepen roots, but this will also increase photosynthesis and vegetation
  !! biomass slightly, due to more access to soil water. so while use of
  !! deeper roots is desirable in the coupled global model, one may decide
  !! to use ctem's default parameterizarion for stand alone simulations, and
  !! set deeproots to .false. 
  deeproots = .false.

  !     ---------------------------------------------------------------
  !
  !     initialization
  !
  rmatctem(:,:,:)=0.0 !il1:il2,icc,ignd
  etmp(:,:,:)    =0.0 !il1:il2,icc,ignd
  rmatc(:,:,:)   =0.0 !il1:il2,ican,ignd
  ailc(:,:)=0.0       !il1:il2,ican
  saic(:,:)=0.0       !il1:il2,ican
  paic(:,:)=0.0       !il1:il2,ican
  slaic(:,:)=0.0      !il1:il2,ican
  zolnc(:,:)=0.0      !il1:il2,ican
  averough(:,:)=0.0   !il1:il2,ican
  alvisc(:,:)=0.0     !il1:il2,ican
  alnirc(:,:)=0.0     !il1:il2,ican
  cmasvegc(:,:)=0.0   !il1:il2,ican
  sfcancmx(:,:)=0.0   !il1:il2,ican
  sla(:)=0.0          !icc
  useb(:,:)=0.0       !il1:il2,icc
  ailcg(:,:)=0.0      !il1:il2,icc
  ailcb(:,:)=0.0      !il1:il2,icc
  veghght(:,:)=0.0    !il1:il2,icc
  lnrghlth(:,:)=0.0   !il1:il2,icc
  a(:,:)=0.0          !il1:il2,icc
  slai(:,:)=0.0       !il1:il2,icc
  sai(:,:)=0.0        !il1:il2,icc
  bmasveg(:,:)=0.0    !il1:il2,icc
  pai(:,:)=0.0        !il1:il2,icc

  !> Generate the sort index for correspondence between CTEM pfts and the
  !>  values in the parameter vectors
  sort = genSortIndex()

  do 80 j = 1,icc
    do 90 i =  il1, il2
      usealpha(i,j)=alpha(sort(j))
90      continue
80    continue
  !
  !>
  !!------ 1. conversion of leaf biomass into leaf area index -------
  !!
  !! find specific leaf area (sla, m2/kg) using leaf life span
  !!
  do j = 1,icc
    sla(j) = 25.0*(lfespany(sort(j))**(-0.50))
    if(specsla(sort(j)).gt.abszero) sla(j)=specsla(sort(j)) 
  end do

  !>
  !! Convert leaf biomass into lai. brown leaves could have less
  !! lai than the green leaves for the same leaf mass. for now we
  !! assume sla of brown leaves is fracbofg times that of green 
  !! leaves. 
  !!
  !! also find stem area index as a function of stem biomass
  !!
  do 150 j = 1,icc
    do 160 i = il1,il2
      if (fcancmx(i,j) > 0.0) then
        ailcg(i,j) = sla(j) * gleafmas(i,j)
        ailcb(i,j) = sla(j) * bleafmas(i,j) * fracbofg
        ! stem area index
        sai(i,j) = 0.55 * (1.0 - exp(-0.175 * stemmass(i,j))) 
        !> Plant area index is sum of green and brown leaf area indices
        !! and stem area index
        pai(i,j) = ailcg(i,j) + ailcb(i,j) + sai(i,j)

        !> Make class see some minimum pai, otherwise it runs into numerical
        !! problems
        pai(i,j) = max(0.3, pai(i,j))

      endif
160 continue
150   continue

  !> Get fcancmx weighted leaf area index for use by class
  !! needle leaf evg + dcd = total needle leaf    
  !! broad leaf evg + dcd cld + dcd dry = total broad leaf    
  !! crop c3 + c4 = total crop
  !! grass c3 + c4 = total grass
  !! also add brown lai. note that although green + brown
  !! lai is to be used by class for energy and water balance
  !! calculations, stomatal conductance estimated by the 
  !! photosynthesis subroutine is only based on green lai.
  !! that is although both green+brown leaves intercept
  !! water and light, only the green portion photosynthesizes.
  !! also lump stem and plant area indices for class' 4 pfts
  !!
  do 200 j = 1, ican
    do 210 m = reindexPFTs(j,1), reindexPFTs(j,2)    
      do 220 i = il1, il2
        sfcancmx(i,j) = sfcancmx(i,j) + fcancmx(i,m)
        ailc(i,j) = ailc(i,j) + (fcancmx(i,m) * (ailcg(i,m) + ailcb(i,m)))
        saic(i,j) = saic(i,j) + (fcancmx(i,m) * sai(i,m))
        paic(i,j) = paic(i,j) + (fcancmx(i,m) * pai(i,m))
        slaic(i,j) = slaic(i,j) + (fcancmx(i,m) * slai(i,m))
220   continue
210 continue
200 continue

  do 230 j = 1, ican
    do 240 i = il1, il2
      if (sfcancmx(i,j) > abszero) then
         ailc(i,j) = ailc(i,j) / sfcancmx(i,j)
         saic(i,j) = saic(i,j) / sfcancmx(i,j)
         paic(i,j) = paic(i,j) / sfcancmx(i,j)
         slaic(i,j) = slaic(i,j) / sfcancmx(i,j)
      else
         ailc(i,j) = 0.0
         saic(i,j) = 0.0
         paic(i,j) = 0.0
         slaic(i,j) = 0.0
      endif

      !>for crops and grasses set the minimum lai to a small number, other
      !!wise class will never run tsolvc and thus phtsyn and ctem will not
      !!be able to grow crops or grasses.

      select case(classpfts(j))
        case('Crops','Grass')
          ailc(i,j) = max(ailc(i,j), 0.1)
        case ('NdlTr' , 'BdlTr', 'BdlSh') 
          ! Do nothing for non-grass/crop 
        case default
          print*,'Unknown CLASS PFT in bio2str ',classpfts(j)
          call XIT('bio2str',-1)                                         
      end select
240 continue
230 continue

  !>------ 2. conversion of stem biomass into roughness length -------
  !!
  !! CLASS uses log of roughness length (zoln) as an input parameter. when 
  !! vegetation grows and dies as per ctem, then zoln is provided by ctem.
  !!
  !! 1. convert stem biomass into vegetation height for trees and crops,
  !! and convert leaf biomass into vegetation height for grass
  !!
  !! 2. convert vegetation height into roughness length & take its log
  !!
  !! 3. lump this for ctem's 9 pfts into class' 4 pfts
  !!

  !> Determine the vegetation height of the peatland shrubs
  !! shrubs in mbbog maximum 0.3 m average 0.18 m (Bubier et al. 2011)
  !! grass and herbs avg 0.3m, min 0.05 max 0.8
  !! low shrub avg 0.82, min 0.10, max 2.00 
  !! tall shrub avg.3.76, min 2.3, max 5.00 (Hopkinson et al. 2005)
  !! The Canadian wetland vegetation classification
  !! graminoids include grass,rush,reed,sedge
  !! forb is all non-graminoids herbaceous plants
  !! shrubs dwarf <0.1m, low (0.1 to 0.5), medium 0.5 to 1.5, tall >1.5 
  !! trees > 5 m  
  do 250 j = 1, ican
    do 260 m = reindexPFTs(j,1), reindexPFTs(j,2)    
      do 270 i = il1, il2
        select case (classpfts(j))
          case ('NdlTr','BdlTr') !Trees
            if (ipeatland(i) .eq. 0) then ! For uplands:
              veghght(i,m)= min(10.0 * stemmass(i,m)**0.385, 45.0)
            else !peatland trees have a different relation than normal. Max height 10 m.
              veghght(i,m) = min(3.0 * stemmass(i,m)**0.385, 10.0) 
            end if
          case ('BdlSh')
            if (ipeatland(i) .eq. 0) then ! For uplands:     
              veghght(i,m) = min(10.0 * stemmass(i,m)**0.385, 45.0) !FLAG SET TO TREES!!!!!!
            else ! peatland shrubs         
              veghght(i,m) = min(1.0, 0.25 * (stemmass(i,m)**0.2))  
            end if
          case ('Crops') ! <Crops
            veghght(i,m) = 1.0 * (stemmass(i,m) + gleafmas(i,m))**0.385
          case ('Grass') ! <Grass
            if (ipeatland(i) .eq. 0) then ! For uplands:     
              veghght(i,m) = 3.5 * (gleafmas(i,m) + fracbofg * bleafmas(i,m))**0.50   
            else ! peatland grasses and sedges                  
              veghght(i,m) = min(1.0, (gleafmas(i,m) + fracbofg * bleafmas(i,m))**0.3) 
            end if
          case default
            print*,'Unknown CLASS PFT in bio2str ',classpfts(j)
            call XIT('bio2str',-2)                                                         
        end select
        lnrghlth(i,m) = log(0.10 * max(veghght(i,m), 0.10))
        
        averough(i,j) = averough(i,j) + (fcancmx(i,m) * lnrghlth(i,m))
        
270   continue
260 continue
250 continue

  do 330 j = 1, ican
    do 340 i = il1, il2
      if(sfcancmx(i,j).gt.abszero)then
        averough(i,j) = averough(i,j) / sfcancmx(i,j)
      else
        averough(i,j) = -4.605
      endif
      zolnc(i,j) = averough(i,j)
340 continue
330 continue

  !>------ 3. estimating fraction of roots in each soil layer for -----
  !!------      ctem's each vegetation type, using root biomass   -----
  !!
  !! Estimate parameter b of variable root profile parameterization
  !!
  do j = 1, icc 
    if(deeproots) then
        b(j) = (abar(sort(j)) - 1.5) * (avertmas(sort(j))**alpha(sort(j)))
    else
        b(j) = abar(sort(j)) * (avertmas(sort(j))**alpha(sort(j)))
    endif
  end do

  !> Use b to estimate 99% rooting depth
  !!
  do 370 j = 1,ican
    do 380 m = reindexPFTs(j,1), reindexPFTs(j,2)    
      do 390 i = il1, il2

        useb(i,m) = b(m)
        usealpha(i,m) = alpha(sort(m))
        rootdpth(i,m) = (4.605 * (rootmass(i,m)**alpha(sort(m)))) / b(m)
        
        !> If estimated rooting depth is greater than the perennially-frozen
        !! soil depth, the permeable soil depth or the maximum rooting depth
        !! then adjust rooting depth and parameter alpha. The soildepth is
        !! the permeable soil depth from the initialization file, zbotw is the
        !! bottom of the soil layer that soil depth lies within.
        !!
        !! Also find "a" (parameter determining root profile). this is
        !! the "a" which depends on time varying root biomass 
        !!
        if(rootdpth(i,m) > min(soildpth(i), maxAnnualActLyr(i), &
                               zbotw(i,ignd), mxrtdpth(sort(m)))) then
                               
          rootdpth(i,m) = min(soildpth(i),maxAnnualActLyr(i),&
                             zbotw(i,ignd),mxrtdpth(sort(m)))
          if (rootdpth(i,m) <= abszero) then
            a(i,m) = 100.0
          else
            a(i,m) = 4.605 / rootdpth(i,m)
          end if
        else
        
          if (rootmass(i,m) <= abszero) then
            a(i,m) = 100.0
          else
            a(i,m) = useb(i,m) / (rootmass(i,m)**usealpha(i,m))
          end if
        end if
390   continue
380 continue
370 continue

  do 400 j = 1,icc
    do 410 i = il1, il2

      kend = 9999  ! initialize with a dummy value

      !> Using parameter "a" we can find fraction of roots in each soil layer just like class
      zroot = rootdpth(i,j)
      totala(i,j) = 1.0 - exp(-a(i,j) * zroot)

      !! if rootdepth is shallower than the bottom of the first layer
      if (zroot <= zbotw(i,1)) then
        rmatctem(i,j,1)=1.0
        do k = 2, ignd
         rmatctem(i,j,k) = 0.0
        end do
        kend = 1
      else
        do k = 2, ignd
          if (zroot <= zbotw(i,k) .and. zroot > zbotw(i,k-1)) then

            !>if rootdepth is shallower than the bottom of current layer and
            !!is deeper than bottom of the previous top layer
            kend = k ! kend = soil layer number in which the roots end  
          end if
        end do

        if (kend .eq. 9999) then
          write(6,2100) i,j,k,kend
2100      format(' at (i) = (',i3,'), pft=',i2,', depth=',i2,' kend &
          & is not assigned. kend  = ',i5)
          call xit('bio2str',-3)
        end if

        etmp(i,j,1) = exp(-a(i,j) * zbotw(i,1))
        rmatctem(i,j,1) = (1.0 - etmp(i,j,1)) / totala(i,j)            
      
        if (kend == 2) then
          !>if rootdepth is shallower than the bottom of 2nd layer
          etmp(i,j,kend) = exp (-a(i,j) * zroot)
          rmatctem(i,j,kend) = (etmp(i,j,kend-1) - etmp(i,j,kend)) &
                               / totala(i,j)
        else if (kend > 2) then
          !>if rootdepth is shallower than the bottom of 3rd layer 
          !!or even the deeper layer (ignd>3)
          do k=2, kend-1
            etmp(i,j,k) = exp(-a(i,j) * zbotw(i,k))            
            rmatctem(i,j,k) = (etmp(i,j,k-1) - etmp(i,j,k)) / totala(i,j)
          end do

          etmp(i,j,kend) = exp(-a(i,j) * zroot)
          rmatctem(i,j,kend) = (etmp(i,j,kend-1) - etmp(i,j,kend)) &
                               / totala(i,j)
        end if   !if kend
      end if    !zroot
410 continue
400 continue


  !> We only allow roots in non-perennially frozen soil layers so first check which layers are
  !! unfrozen and then adjust the distribution appropriately. For defining which
  !! layers are frozen, we use the active layer depth.
  rmatctem = unfrozenRoots(il1,il2,ilg,maxAnnualActLyr,zbotw,rmatctem)

  !> Make sure all fractions (of roots in each layer) add to one.
  do 411 j = 1, icc
    do 412 i = il1, il2
      if (fcancmx(i,j) > 0.0) then
        rmat_sum = 0.0
        do k = 1, ignd
          rmat_sum = rmat_sum + rmatctem(i,j,k)
        end do
        if( abs(rmat_sum-1.0).gt.1e-10) then
          write(6,2300) i,j,rmat_sum
2300      format(' at (i) = (',i3,'), pft=',i2,' fractions of roots&
            not adding to one. sum  = ',f12.7)
          call xit('bio2str',-4)
        end if
      end if
412 continue
411 continue

  !> Lump rmatctem(i,icc,ignd)  into rmatc(i,ican,ignd) for use by CLASS
  do j = 1, ican
    do m = reindexPFTs(j,1), reindexPFTs(j,2)    
      do i = il1, il2
        do k = 1, ignd
          rmatc(i,j,k) = rmatc(i,j,k) + (fcancmx(i,m) * rmatctem(i,m,k))  
        end do
      end do
    end do      
  end do        

  do 450 j = 1, ican
    do 460 i = il1, il2
      if (sfcancmx(i,j) > abszero) then
        do k = 1, ignd
          rmatc(i,j,k) = rmatc(i,j,k) / sfcancmx(i,j)
        end do
      else
        rmatc(i,j,1) = 1.0
        do k = 2, ignd
          rmatc(i,j,k) = 0.0
        end do
      end if
460 continue
450 continue

  !>-------------------  4. calculate storage lai  --------------------
  do 500 j = 1, icc
    do 510 i = il1, il2
      if (fcancmx(i,j) > 0.0) then
        slai(i,j) = ((stemmass(i,j) + rootmass(i,j)) &
                    / eta(sort(j)))**(1. / kappa(sort(j)))
        slai(i,j) = (prcnslai(sort(j)) / 100.0) * sla(j) * slai(i,j)

        !>need a minimum slai to be able to grow from scratch. consider this as model seeds.
        slai(i,j) = max(slai(i,j), minslai(sort(j)))
      end if
510 continue
500 continue
  
  !>--- 5. calculate total vegetation biomass for each ctem pft, and --
  !!---------------- canopy mass for each class pft ------------------
  do 550 j = 1, icc
    do 560 i = il1, il2
      if (fcancmx(i,j) > 0.0) then
        bmasveg(i,j) = gleafmas(i,j) + stemmass(i,j) + rootmass(i,j)
      end if
560 continue
550 continue

  !!Since class uses canopy mass and not total vegetation biomass as an
  !!input, we find canopy mass as a sum of stem and leaf mass, for each
  !!class pft, i.e. only above ground biomass. 
  do 600 j = 1, ican
    do 610 m = reindexPFTs(j,1), reindexPFTs(j,2)    
      do 620 i = il1, il2
        cmasvegc(i,j)= cmasvegc(i,j) + (fcancmx(i,m) &
                       * (bleafmas(i,m) + gleafmas(i,m) + stemmass(i,m)))  
620   continue
610 continue
600 continue

  do 630 j = 1, ican
    do 640 i = il1, il2
      if (sfcancmx(i,j) > abszero) then
        cmasvegc(i,j) = cmasvegc(i,j) / sfcancmx(i,j)
        cmasvegc(i,j) = cmasvegc(i,j) * (1.0 / 0.50) !assuming biomass is 50% C
      else
        cmasvegc(i,j) = 0.0
      end if

      !> If there is no vegetation canopy mass will be abszero. this should 
      !!essentially mean more bare ground, but since we are not changing
      !!fractional coverages at present, we pass a minimum canopy mass
      !!to class so that it doesn't run into numerical problems.
      !          cmasvegc(i,j)=max(cmasvegc(i,j),3.0)    !YW April 14, 2015 ! FLAG JM - is this okay? Nov 2016.
      cmasvegc(i,j) = max(cmasvegc(i,j), 0.1)     
640 continue
630 continue

  !>--- 6. calculate albedo for class' 4 pfts based on specified ----
  !!------ albedos of ctem 9 pfts and their fractional coveraes -----
  do 700 j = 1, ican
    do 710 m = reindexPFTs(j,1), reindexPFTs(j,2)    
      do 720 i = il1, il2
        alvisc(i,j) = alvisc(i,j) + (fcancmx(i,m) * albvis(sort(m)))  
        alnirc(i,j) = alnirc(i,j) + (fcancmx(i,m) * albnir(sort(m)))  
720   continue
710 continue
700 continue

  do 730 j = 1, ican
    do 740 i = il1, il2
      if (sfcancmx(i,j) > abszero) then
        alvisc(i,j) = (alvisc(i,j) / sfcancmx(i,j)) / 100.0
        alnirc(i,j) = (alnirc(i,j) / sfcancmx(i,j)) / 100.0
      else
        alvisc(i,j) = 0.0
        alnirc(i,j) = 0.0
      end if
740 continue
730 continue

  return
  
  end subroutine allometry
!!@}
! ---------------------------------------------------------------------------------------------------

!>\ingroup allometry_unfrozenRoots
!!@{ We only allow roots in non-perennially frozen soil layers so first check which layers are
!! unfrozen and then correct the root distribution appropriately. For defining which
!! layers are frozen, we use the active layer depth.
!!@author J. Melton
function unfrozenRoots(il1,il2,ilg,maxAnnualActLyr,zbotw,rmatctem)

  use classic_params, only : ignd,icc
  
  implicit none

  integer, intent(in) :: il1             !< il1=1
  integer, intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
  integer, intent(in) :: ilg             !< Number of grid cells/tiles in latitude circle
  real, intent(in)    :: maxAnnualActLyr(:)!< Active layer depth maximum over the e-folding period specified by parameter eftime (m).
  real, intent(in)    :: zbotw(:,:)      !< Bottom of soil layers (m)
  real, intent(in)    :: rmatctem(:,:,:) !<fraction of roots for each of ctem's pfts in each soil layer

  real :: unfrozenRoots(ilg,icc,ignd)             !< root distribution only over unfrozen layers  
  
  integer :: botlyr                       !< bottom layer of the unfrozen soil column
  real :: frznrtlit                       !< fraction of root distribution in frozen layers
  integer k,j,i

  !! We only add to non-perennially frozen soil layers so first check which layers are
  !! unfrozen and then do the allotment appropriately. For defining which
  !! layers are frozen, we use the active layer depth.
  unfrozenRoots = 0.
  do i = il1, il2
    ! Find the bottom of the unfrozen soil column:
    botlyr = 1 ! we assume if the first layer is frozen that it still can
                 ! accept the root litter. So initialize to 1.
    do k = 1,ignd
        if(maxAnnualActLyr(i) < zbotw(i,k)) exit
          botlyr = k
    end do

    do j = 1, icc      
      if (botlyr == ignd) then !if the botlyr is the bottom of the soil column then just set to original
                               ! and move on.
           unfrozenRoots(i,j,:) = rmatctem(i,j,:)
      else ! there is some frozen soil so adjust how the root litter is distibuted
        frznrtlit = sum(rmatctem(i,j,botlyr+1:ignd)) !determine how much of the distribution is in the frozen layers
        do k = 1, botlyr
            unfrozenRoots(i,j,k) = rmatctem(i,j,k) + rmatctem(i,j,k) / (1. - frznrtlit) * frznrtlit
        end do
      end if
    end do 
  end do 
  
end function unfrozenRoots
!>@}
! -----------------------------------------------------------------------------------------------
!>\namespace allometry
!!

!! Converts biomass to structural attributes
!!@author V. Arora, J. Melton, Y. Peng 
!!
!!The time-varying biomass in the leaves (\f$C_L\f$), stem (\f$C_S\f$) and root
!! (\f$C_R\f$) components is used to calculate the structural attributes of vegetation
!! for the energy and water balance calculations by CLASS.
!!
!!Leaf biomass is converted to LAI using specific leaf area
!! (\f${SLA}\f$, \f$m^2\,(kg\,C)^{-1}\f$), which itself is assumed
!! to be a function of leaf lifespan (\f$\tau_L\f$; see also classic_params.f90)
!!
!!\f[ \label{sla} SLA= \gamma_L\tau_L^{-0.5}\\ LAI = C_LSLA\nonumber \f]
!!
!!where \f$\gamma_L\f$ is a constant with value equal to \f$25\,m^2\,(kg\,C)^{-1}\,
!!yr^{0.5}\f$. The vegetation height (\f$H\f$; \f$m\f$) is calculated for tree, crop
!! and grass PFTs as
!!
!!\f$ H = \min (10.0C_S^{0.385},45) \f$ for trees
!!
!!\f$ H = (C_S + C_L)^{0.385} \f$ for crops
!!
!!\f$ H = 3.5 (C_{L,g} + 0.55C_{L,b})^{0.5} \f$ for  grasses
!!
!!where \f$C_{L,g}\f$ is the green leaf biomass and \f$C_{L,b}\f$ is the brown
!! leaf biomass that is scaled by 0.55 to reduce its contribution to the plant
!! height. CTEM explicitly tracks brown leaf mass for grass PFTs. The turnover of
!! green grass leaves, due to normal aging or stress from drought and/or cold, does
!! not contribute to litter pool directly as the leaves first turn brown. The brown
!! leaves themselves turnover to litter relatively rapidly \f$(\tau_{L,b} = 0.1\,\tau_L\f$).
!!
!! In peatlands the vegetation height is calculated as (Wu et al. 2016) \cite Wu2016-zt
!!
!!\f$ H = \min(10.0,3.0C_S^{0.385}) \f$ for trees
!!
!!\f$ H = \min(1.0, 0.25(C_S^{0.2})) \f$ for shrubs
!!
!!\f$ H = = min(1.0,(C_{L,g}+0.55C_{L,b})^{0.3}) \f$ for grasses and sedges
!!
!!CTEM dynamically simulates root distribution and depth in soil following (Arora and Boer, 2003) \cite Arora2003838. The root distribution takes an exponential form and roots grow and deepen with increasing root biomass. The cumulative root fraction at depth \f$z\f$ is given by
!!
!!\f$ f_R(z) = 1 - \exp(-\iota z) \f$
!!
!!Rooting depth (\f$d_R\f$; \f$m\f$), which is defined to be the depth 
!!containing \f$99\,{\%}\f$ of the root mass, is found by setting \f$z\f$ equal 
!!to \f$d_R\f$ and \f$f_R = 0.99\f$, which yields
!!
!!\f[ \label{rootterm1} d_R = \frac{-\ln(1-f_R)}{\iota} = \frac{-\ln(1 - 0.99)}{\iota} = \frac{4.605}{\iota}. \f]
!!
!!The parameter \f$\iota\f$ that describes the exponential root distribution is calculated as
!!\f[ \label{iota} \iota = \overline{\iota} \left(\frac{\overline{C_R}}{C_R} \right)^{0.8}, \f]
!!
!!where \f$\overline{\iota}\f$ represents the PFT-specific mean root distribution 
!! profile parameter and \f$\overline{C_R}\f$ the average root biomass derived 
!! from Jackson et al. (1996) \cite Jackson1996-va (see also classic_params.f90).
!! Equation for \f$\iota\f$ above yields a lower (higher) value of \f$\iota\f$ than
!! \f$\overline{\iota}\f$ when root biomass \f$C_R\f$ is higher (lower) than the 
!!PFT-specific mean root biomass \f$\overline{C_R}\f$, resulting in a deeper
!! (shallower) root profile than the mean root profile.
!!
!!The rooting depth \f$d_R\f$ is checked to ensure it does not exceed the soil 
!!depth or extend into perennially frozen soil. If so, \f$d_R\f$ is set to the soil depth or mean
!!annual maximum active layer depth (based on e-folding time of 5 years), whichever is shallower, 
!! and \f$\iota\f$ is recalculated
!! as \f$\iota = 4.605/d_R\f$ (see Eq. \ref{rootterm1} for derivation of 4.605 term).
!! The new value of \f$\iota\f$ is used to determine the root distribution profile 
!!adjusted to the shallower depth. Finally, the root distribution profile is used 
!!to calculate fraction of roots in each of the model's soil layers.
!!
!>\file
end module applyAllometry