      subroutine ch4wetland (hetrores, il1, il2, ilg, ta, wetfrac,
     1                             ig, npp, tbar, thliqg, currlat,
     2                          sand,  wetfrac_s,
c
c    -------------- inputs above this line, outputs below -------------
     2                        ch4wet1,    ch4wet2,    wetfdyn,  
     3                        ch4dyn1,    ch4dyn2)
c
c               canadian terrestrial ecosystem model (ctem) v1.0
c           	methane flux from wetlands subroutine 
c
c     09  june. 2010 - this subroutine calculates methane flux
c     Brian Amiro      from wetlands of a grid cell (i.e. land only)

c     inputs 
c
c     hetrores  - heterotrophic respiration from main ctem program
c				calculated as sum of litres + socres
c     ilg       - no. of grid cells in latitude circle
c     il1,il2   - il1=1, il2=ilg
c     ig        - no. of soil layers (3)
c     ta        - air temperature, k
c     wetfrac   - prescribed fraction of wetlands in a grid cell
c     npp       - grid-averaged npp from ctem (u-mol co2/m2.s)
c     tbar      - temperature of soil layers
c     thliqg    - liquid soil moisture content (fraction)
c     currlat   - centre latitude of grid cells in degrees
c     sand      - percentage sand in soil layers
c     wetfrac_s - prescribed fraction of wetlands based on slope only 
c
c     outputs
c
c     ch4wet1   - methane flux from wetlands calculated using hetrores 
c                 in umol ch4-c/m2.s
c     ch4wet2   - methane flux from wetlands calculated using npp
c                 in umol ch4-c/m2.s
c     wetfdyn   - dynamic gridcell wetland fraction determined using 
c                 slope and soil moisture
c     
c     ch4dyn1   - methane flux from wetlands calculated using hetrores
c                 and wetfdyn, in umol ch4-c/m2.s
c     ch4dyn2   - methane flux from wetlands calculated using npp 
c                 and wetfdyn, in umol ch4-c/m2.s
c
      use ctem_params,        only : wtdryres, ratioch4, factor2

      implicit none
c
      integer ilg, il1, il2, i, ig
c
      real      ch4wet1(ilg),  wetresp(ilg), ta(ilg), wetfrac (ilg),
     1          ch4dyn1(ilg),  ch4dyn2(ilg)
      real     hetrores(ilg),  ch4wet2(ilg),        wetfrac_s(ilg)
      real          npp(ilg),      tbar(ilg,ig),    thliqg(ilg, ig),
     1          sand(ilg,ig),      wetfdyn(ilg),       currlat(ilg),
     2              porosity,      soil_wetness
c
c     ---------------------------------------------------------------
c     Constants and parameters are located in ctem_params.f90
c     -----------------------------------------------------------------
c
c     initialize required arrays to zero
c
      do 110 i = il1, il2
        wetresp(i)=0.0     ! heterotrophic wetland respiration
        ch4wet1(i)=0.0      ! methane flux from wetlands
        ch4wet2(i)=0.0
        ch4dyn1(i)=0.0
        ch4dyn2(i)=0.0
110   continue
c
c
c     initialization ends  
c	--------------------------------------------------
c
c	estimate the methane flux from wetlands for each grid cell
c	scaling by the wetland fraction in a grid cell
c 	and set the methane flux to zero when screen temperature (ta) is 
c	below or at freezing
c	this is consistent with recent flux measurements by the university
c	of manitoba at churchill, manitoba
c
c     debugging, apr 15/2011
c        call xit('ch4wtlnd',-1)
c
	do 210 i = il1, il2 
           wetresp(i)=hetrores(i)*wtdryres*wetfrac(i)
	   ch4wet1(i)=ratioch4*wetresp(i)
           ch4wet2(i)=factor2*wetfrac(i)*max(0.0,npp(i))
     1                *(2**((tbar(i,1)-273.2)/10.0)) 
           if (ta(i).lt.273.2) then
             ch4wet1(i)=0.0
             ch4wet2(i)=0.0
           endif
c
           porosity=(-0.126*sand(i,1)+48.9)/100.0 ! top soil layer porosity
           soil_wetness=(thliqg(i,1)/porosity)
           soil_wetness=max(0.0,min(soil_wetness,1.0))
c
c          if soil wetness meets a latitude specific threshold then the slope
c          based wetland fraction is wet and is an actual wetland else not
c
           wetfdyn(i)=0.0  ! initialize dynamic wetland fraction to zero
c
           if (currlat(i).ge.50.0) then ! all area north of 35 n
c            if(soil_wetness.gt.0.635)then ! cku
c            if(soil_wetness.gt.0.55)then ! ckw 3031-3040
             if(soil_wetness.gt.0.60)then ! ckw 3041-3050
               wetfdyn(i)=wetfrac_s(i)
             endif
           elseif (currlat(i).lt.50.0.and.currlat(i).ge.-10.0) then ! between 10 s and 35 n
c            if(soil_wetness.gt.0.84)then   ! cku
c            if(soil_wetness.gt.0.75)then   ! ckw 3031-3040
             if(soil_wetness.gt.0.80)then   ! ckw 3041-3050
               wetfdyn(i)=wetfrac_s(i)
             endif
           else ! everything else below 10 s
c            if(soil_wetness.gt.0.75)then !cku 
c            if(soil_wetness.gt.0.68)then ! ckw 3031-3040
             if(soil_wetness.gt.0.70)then ! ckw 3041-3050
               wetfdyn(i)=wetfrac_s(i)
             endif
           endif           

c new dynamic calculation
c same as ch4wet1 & 2, but wetfrac replaced by wetfdyn
           wetresp(i)=hetrores(i)*wtdryres*wetfdyn(i)
           ch4dyn1(i)=ratioch4*wetresp(i)
           ch4dyn2(i)=factor2*wetfdyn(i)*max(0.0,npp(i))
     1                *(2**((tbar(i,1)-273.2)/10.0))
           if (ta(i).lt.273.2) then
             ch4dyn1(i)=0.0
             ch4dyn2(i)=0.0
           endif

210	continue

        return
        end
 




