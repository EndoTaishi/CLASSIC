
       subroutine  balcar (gleafmas, stemmass, rootmass, bleafmas,      
     1                     litrmass, soilcmas, ntchlveg, ntchsveg,
     2                     ntchrveg, tltrleaf, tltrstem, tltrroot,
     3                     glcaemls, blcaemls, stcaemls, rtcaemls,
     4                     ltrcemls, ltresveg, scresveg, humtrsvg,
     5                     pglfmass, pblfmass, pstemass, protmass,
     6                     plitmass, psocmass, vgbiomas, repro_cost,
     7                     pvgbioms, gavgltms, pgavltms, gavgscms,
     8                     pgavscms, galtcels, repro_cost_g,
     9                          npp,  autores, hetrores,      gpp,
     a                          nep,   litres,   socres, dstcemls,
     b                          nbp, litrfall, humiftrs,
     c                          il1,      il2
c    --------------moss variables--------------------------------------\
C
     1              ,ipeatland, Cmossmas, pCmossmas
     2              ,nppmosstep, litrfallms, litrmassms,plitrmassms,
     3              ltrestepms,humicmstep,socrestep,hutrstep)    
c    -----------YW March 20, 2015 -------------------------------------/   

c     -----------------------------------------------------------------      
c
c               Canadian Terrestrial Ecosystem Model (CTEM)
C                          Carbon Balance Subroutine
c
c     22  Nov 2012  - calling this version 1.1 since a fair bit of ctem
c     V. Arora        subroutines were changed for compatibility with class
c                     version 3.6 including the capability to run ctem in
c                     mosaic/tile version along with class.
c
c     24  Sep 2012  - add in checks to prevent calculation of non-present
c     J. Melton       pfts
c
c     27  May 2003  - this subroutine checks if the various c fluxes
c     V. Arora        between the different pools balance properly to
c                     make sure that conservation of mass is achieved 
c                     with in a specified tolerance.
c     inputs 
c
c                 unless mentioned all pools are in kg c/m2
c
c                 pools (after being updated)
c
c     stemmass  - stem mass for each of the 9 ctem pfts
c     rootmass  - root mass for each of the 9 ctem pfts
c     gleafmas  - green leaf mass for each of the 9 ctem pfts
c     bleafmas  - brown leaf mass for each of the 9 ctem pfts
c     litrmass  - litter mass over the 9 pfts and the bare fraction
c                 of the grid cell
c     soilcmas  - soil carbon mass over the 9 pfts and the bare fraction
c                 of the grid cell
c     repro_cost - amount of C transferred to litter due to reproductive tissues
c
c                 grid averaged pools
c     vgbiomas  - vegetation biomass
c     gavgltms  - litter mass
c     gavgscms  - soil carbon mass
c
c                 pools (before being updated)

c                 variable explanation same as above.
c     pglfmass  - previous green leaf mass
c     pblfmass  - previous brown leaf mass
c     pstemass  - previous stem mass
c     protmass  - previous root mass
c     plitmass  - previous litter mass
c     psocmass  - previous soil c mass
c       
c                 grid average pools
c     pvgbioms  - previous vegetation biomass
c     pgavltms  - previous litter mass
c     pgavscms  - previous soil c mass
c                 
c                 unless mentioned all fluxes are in units of 
c                 u-mol co2/m2.sec
c
c                 fluxes for each pft
c
c     ntchlveg  - net change in leaf biomass
c     ntchsveg  - net change in stem biomass
c     ntchrveg  - net change in root biomass
c                 the net change is the difference between allocation
c                 and autotrophic respiratory fluxes 
c
c     tltrleaf  - total leaf litter falling rate
c     tltrstem  - total stem litter falling rate
c     tltrroot  - total root litter falling rate
c 
c                 carbon emission losses mainly due to fire 
c     glcaemls  - green leaf carbon emission losses 
c     blcaemls  - brown leaf carbon emission losses 
c     stcaemls  - stem carbon emission losses 
c     rtcaemls  - root carbon emission losses 
c     ltrcemls  - litter carbon emission losses
c
c     ltresveg  - litter respiration for each pft + bare fraction 
c     scresveg  - soil c respiration for each pft + bare fraction
c     humtrsvg  - humification for each pft + bare fraction
c
c                 grid averaged fluxes
c
c     npp       - net primary productivity 
c     autores   - autotrophic respiration
c     hetrores  - heterotrophic respiration
c     gpp       - gross primary productivity
c     nep       - net primary productivity     
c     litres    - litter respiration
c     socres    - soil carbon respiration
c     dstcemls  - carbon emission losses due to disturbance, mainly fire
c     galtcels  - carbon emission losses from litter
c     expnbaln  - amount of c related to spatial expansion
c     repro_cost_g - amount of C used to generate reproductive tissues
c     nbp       - net biome productivity
c     litrfall  - combined (leaves, stem, and root) total litter fall rate   
c     humiftrs  - humification
c
c                 other variables
c 
c     deltat    - ctem's time step
c     icc       - no. of ctem plant function types, currently 8
c     ilg       - no. of grid cells in latitude circle
c     il1,il2   - il1=1, il2=ilg
c     fcancmx   - max. fractional coverage of ctem's 9 pfts, but this can be
c                modified by land-use change, and competition between pfts

      use ctem_params,        only : tolrance, icc, ilg, deltat
c
      implicit none
c
      integer il1, il2, i, j, k
c
      real stemmass(ilg,icc),   rootmass(ilg,icc),   gleafmas(ilg,icc),    
     1     bleafmas(ilg,icc), litrmass(ilg,icc+1), soilcmas(ilg,icc+1), 
     2     ntchlveg(ilg,icc),   ntchsveg(ilg,icc),   ntchrveg(ilg,icc),
     3     tltrleaf(ilg,icc),   tltrstem(ilg,icc),   tltrroot(ilg,icc),
     4     glcaemls(ilg,icc),   blcaemls(ilg,icc),   stcaemls(ilg,icc),
     5     rtcaemls(ilg,icc),   ltrcemls(ilg,icc), ltresveg(ilg,icc+1),
     6   scresveg(ilg,icc+1), humtrsvg(ilg,icc+1),   pglfmass(ilg,icc),
     7     pblfmass(ilg,icc),   pstemass(ilg,icc),   protmass(ilg,icc),
     8   plitmass(ilg,icc+1), psocmass(ilg,icc+1),            npp(ilg),
     9         vgbiomas(ilg),       pvgbioms(ilg),       gavgltms(ilg),
     a         pgavltms(ilg),       gavgscms(ilg),       pgavscms(ilg),
     9          autores(ilg),       hetrores(ilg),            gpp(ilg),
     a              nep(ilg),         litres(ilg),         socres(ilg),
     b         dstcemls(ilg),            nbp(ilg),       litrfall(ilg),
     c         humiftrs(ilg),       repro_cost(ilg,icc),          
!     d         galtcels(ilg),       expnbaln(ilg),   
     d         galtcels(ilg),    repro_cost_g(ilg)
c
      real             diff1,               diff2
c
c    ------peatland variables ----------------------------------------\
      integer  ipeatland (ilg)
      real     cmossmas(ilg),      pcmossmas(ilg),
     1         nppmosstep(ilg),    litrfallms(ilg),    
     2         litrmassms(ilg),    plitrmassms(ilg),
     3         ltrestepms(ilg),    humicmstep(ilg),
     4         socrestep(ilg),     hutrstep(ilg,icc+1)
c     ------YW March 27, 2015 -----------------------------------------/  
c
      if(icc.ne.12)                            call xit('balcar',-1)  !YW for peatland
c
c     to check c budget we go through each pool for each vegetation
c     type.
c
c     green and brown leaves
c
      do 100 j = 1, icc
        do 110 i = il1, il2
          diff1=(gleafmas(i,j)+bleafmas(i,j)- pglfmass(i,j)-
     &     pblfmass(i,j))
          diff2=(ntchlveg(i,j)- tltrleaf(i,j)- glcaemls(i,j)- 
     &     blcaemls(i,j))*(deltat/963.62)
          if((abs(diff1-diff2)).gt.tolrance)then
            write(6,2000)i,j,abs(diff1-diff2),tolrance
2000        format('at (i)= (',i3,'), pft=',i2,', ',f12.6,' is greater    
     & than our tolerance of ',f12.6,' for leaves')  
            call xit('balcar',-2)
          endif
c         endif
110     continue
100   continue     
c
c     stem
c
      do 150 j = 1, icc
        do 160 i = il1, il2
          diff1=stemmass(i,j) - pstemass(i,j)
          diff2=(ntchsveg(i,j)- tltrstem(i,j)- 
     &     stcaemls(i,j))*(deltat/963.62)
          if((abs(diff1-diff2)).gt.tolrance)then
            write(6,2001)i,j,abs(diff1-diff2),tolrance
2001        format('at (i)= (',i3,'), pft=',i2,', ',f12.6,' is greater
     & than our tolerance of ',f12.6,' for stem')
            call xit('balcar',-3)
          endif
c         endif
160     continue
150   continue    
c
c     root
c
      do 200 j = 1, icc
        do 210 i = il1, il2
          diff1=rootmass(i,j) - protmass(i,j)
          diff2=(ntchrveg(i,j)- tltrroot(i,j)- 
     &     rtcaemls(i,j))*(deltat/963.62)
          if((abs(diff1-diff2)).gt.tolrance)then
            write(6,2002)i,j,abs(diff1-diff2),tolrance
2002        format('at (i)= (',i3,'), pft=',i2,', ',f12.6,' is greater
     & than our tolerance of ',f12.6,' for root')
            call xit('balcar',-4)
          endif
c         endif
210     continue
200   continue    
c
c     litter over all pfts
c
      do 250 j = 1, icc
        do 260 i = il1, il2
          diff1=litrmass(i,j) - plitmass(i,j)
          diff2=( tltrleaf(i,j)+tltrstem(i,j)+tltrroot(i,j)-
     &      ltresveg(i,j)-humtrsvg(i,j)-ltrcemls(i,j)
     &      + repro_cost(i,j))*(deltat/963.62)   
          if((abs(diff1-diff2)).gt.tolrance)then
            write(6,2003)i,j,abs(diff1-diff2),tolrance
2003        format('at (i)= (',i3,'), pft=',i2,', ',f12.6,' is greater
     & than our tolerance of ',f12.6,' for litter')
            call xit('balcar',-5)
          endif
c         endif
260     continue
250   continue    
c
c     litter over the bare fraction
c
        do 280 i = il1, il2
           if (ipeatland(i)==0)         then      !YW May 11, 2015 
          diff1=litrmass(i,icc+1) - plitmass(i,icc+1)
          diff2=( -ltresveg(i,icc+1)-humtrsvg(i,icc+1))*
     &          ( deltat/963.62 )  
          if((abs(diff1-diff2)).gt.tolrance)then
            write(6,2003)i,icc+1,abs(diff1-diff2),tolrance
            call xit('balcar',-6)
          endif
           endif
280     continue
c
c
c     soil carbon over the bare fraction
c
      do 300 j = 1, icc+1
        do 310 i = il1, il2
         if (ipeatland(i)==0)          then
          diff1=soilcmas(i,j) - psocmass(i,j)
          diff2=( humtrsvg(i,j)-scresveg(i,j) )*(deltat/963.62)  
          if((abs(diff1-diff2)).gt.tolrance)then
            write(6,3001)'soilCmas(',i,')=',soilcmas(i,j)
            write(6,3001)'psocmass(',i,')=',psocmass(i,j)
            write(6,3001)'humtr(',i,')=',humtrsvg(i,j)*(deltat/963.62)
            write(6,3001)'scres(',i,')=',scresveg(i,j)*(deltat/963.62)
            write(6,2004)i,j,abs(diff1-diff2),tolrance
2004        format('at (i)= (',i3,'), pft=',i2,', ',f12.6,' is greater
     & than our tolerance of ',f12.6,' for soil c')
            call xit('balcar',-7)
          endif
         endif
310     continue
300   continue   
c
c    --------------------add moss C balance----------------------------\

       do 400 i = il1, il2
           if (ipeatland(i).gt. 0)      then
               diff1 = Cmossmas(i)- pCmossmas(i)
               diff2 = nppmosstep(i) - litrfallms(i)
               if((abs(diff1-diff2)).gt.tolrance)then
                   write(6,3001)'Cmossmas(',i,')=',Cmossmas(i)
                   write(6,3001)'pCmossmas(',i,')=',pCmossmas(i)
                   write(6,3001)'nppmosstep(',i,')=',nppmosstep(i)
                   write(6,3001)' litrfallms(',i,')=',litrfallms(i)
                   write(6,2008)i,abs(diff1-diff2),tolrance
2008           format('at (i)= (',i3,'),',f12.6,' is greater'
     1         'than our tolerance of ',f12.6,' for moss carbon')
                    call xit('balcar',-11)
               endif
               
c    -------------------add moss litter pool C balance-----------------
               
               diff1 = litrmassms(i)- plitrmassms(i)
               diff2 = litrfallms(i)-ltrestepms(i)-humicmstep(i)
               if((abs(diff1-diff2)).gt.tolrance)then
                   write(6,3001)'litrmassms(',i,')=',litrmassms(i)
                   write(6,3001)'plitrmassms(',i,')=',plitrmassms(i)
                   write(6,3001)'litrfallms(',i,')=',litrfallms(i)
                   write(6,3001)' ltrestepms(',i,')=',ltrestepms(i)
                   write(6,3001)' humicmstep(',i,')=',humicmstep(i)
                   write(6,2009)i,abs(diff1-diff2),tolrance
2009               format('at (i)= (',i3,'),',f12.6,' is greater
     1             than our tolerance of ',f12.6,' for moss litter')
                   call xit('balcar',-12)
               endif
          endif
400     continue
c    -------------------YW March 27, 2015------------------------------/
c
c     grid averaged fluxes must also balance
c
c     vegetation biomass
c
      do 350 i = il1, il2
        diff1=vgbiomas(i)-pvgbioms(i)
        diff2=(gpp(i)-autores(i)-litrfall(i)-
     &   dstcemls(i)-repro_cost_g(i))*(deltat/963.62)
        if((abs(diff1-diff2)).gt.(tolrance+0.00003))   then !YW the difference for moss litter 
                                                            !and biomass can go a bit over the tolerance
          write(6,3001)'vgbiomas(',i,')=',vgbiomas(i)
          write(6,3001)'pvgbioms(',i,')=',pvgbioms(i)
          write(6,3001)'     gpp(',i,')=',gpp(i)
          write(6,3001)' autores(',i,')=',autores(i)
          write(6,3001)'litrfall(',i,')=',litrfall(i)
          write(6,3001)'dstcemls(',i,')=',dstcemls(i)
          write(6,3001)'repro_cost_g(',i,')=',repro_cost_g(i)
3001      format(a9,i2,a2,f14.9) 
          write(6,2005)i,abs(diff1-diff2),tolrance
2005      format('at (i)= (',i3,'),',f12.6,' is greater
     & than our tolerance of ',f12.6,' for vegetation biomass')
          write(90,*)    abs(diff1-diff2),tolrance
          call xit('balcar',-8)
        endif
350   continue
c
c     litter
c
      do 380 i = il1, il2
        diff1=gavgltms(i)-pgavltms(i)
        diff2=(litrfall(i)-litres(i)-humiftrs(i)-galtcels(i)
     &   +repro_cost_g(i))*
     &   (deltat/963.62)
        if((abs(diff1-diff2)).gt.(tolrance+0.00003))then    !YW the difference for moss litter 
                                                            !and biomass can go a bit over the tolerance 
          write(6,3001)'pgavltms(',i,')=',pgavltms(i)
          write(6,3001)'gavgltms(',i,')=',gavgltms(i)
          write(6,3001)'litrfall(',i,')=',litrfall(i)
          write(6,3001)'  litres(',i,')=',litres(i)
          write(6,3001)'humiftrs(',i,')=',humiftrs(i)
          write(6,3001)'galtcels(',i,')=',galtcels(i)
          write(6,2006)i,abs(diff1-diff2),tolrance
          write(*,*)i,abs(diff1-diff2),tolrance
2006      format('at (i)= (',i3,'),',f12.6,' is greater
     & than our tolerance of ',f12.6,' for litter mass')
          call xit('balcar',-9)
        endif
380   continue
c
c     soil carbon
c
      do 390 i = il1, il2
        diff1=gavgscms(i)-pgavscms(i)
        diff2=(humiftrs(i)-socres(i))*(deltat/963.62)
        if((abs(diff1-diff2)).gt.tolrance)then
          write(6,3001)'pgavscms(',i,')=',pgavscms(i)
          write(6,3001)'gavgscms(',i,')=',gavgscms(i)
          write(6,3001)'humiftrs(',i,')=',humiftrs(i)
          write(6,3001)'  socres(',i,')=',socres(i)
          write(6,2007)i,abs(diff1-diff2),tolrance
2007      format('at (i)= (',i3,'),',f12.6,' is greater
     & than our tolerance of ',f12.6,' for soil c mass')
          call xit('balcar',-10)
        endif
390   continue
c
      return
      end
