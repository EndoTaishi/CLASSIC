module landuse_change

! Central module for all land use change operations

! J. Melton. Jan 11, 2013

implicit none

real, parameter :: seed = 0.001 ! seed pft fraction, same as in competition
                                ! in mosaic mode, all tiles are given this
                                ! as a minimum

! subroutines contained in this module:
public  :: initialize_luc
public  :: readin_luc
public  :: luc
private :: adjust_luc_fracs

contains

!-------------------------------------------------------------------------------------------------------------
subroutine initialize_luc(nlat,nmos,iyear,lucdat,nmtest,nltest,&
                          mosaic,icc,ican,icp1,nol2pfts,cyclemet,   &
                          cylucyr,lucyr,fcanrow,farerow,nfcancmxrow,  &
                          pfcancmxrow,fcancmxrow,reach_eof)

!           Canadian Terrestrial Ecosystem Model (CTEM) 
!                    LUC Initial Read-In Subroutine 
!
!     9  Jan. 2013  - this subroutine takes in luc information from
!     J. Melton       a luc file and adapts them for runclassctem
!                     it is run once to set up the luc info before the 
!                     model timestepping begins.
!		      

implicit none

! arguments:

! inputs
integer, intent(in) :: nmos
integer, intent(in) :: nlat
integer, intent(in) :: iyear
character(80), intent(in) :: lucdat
logical, intent(in) :: mosaic
integer, intent(in) :: nmtest
integer, intent(in) :: nltest
integer, intent(in) :: icc
integer, intent(in) :: ican
integer, dimension(ican), intent(in) :: nol2pfts
integer, intent(in) :: icp1
logical, intent(in) :: cyclemet
integer, intent(in) :: cylucyr 


! updates
real, dimension(nlat,nmos,icp1), intent(inout) :: fcanrow
real, dimension(nlat,nmos),      intent(inout) :: farerow
logical, intent(inout) :: reach_eof

! outputs
real, dimension(nlat,nmos,icc), intent(out) :: pfcancmxrow
real, dimension(nlat,nmos,icc), intent(out) :: fcancmxrow
real, dimension(nlat,nmos,icc), intent(out) :: nfcancmxrow
integer, intent(out) :: lucyr

! local variables:
real, dimension(icc) :: temparray
real :: temp
integer :: j,m,i,n
integer :: k2,k1,strlen

!-------------------------

!       reset the composite fcanrow as it is appended on later in a loop
        if (.not. mosaic) fcanrow = 0.0

! it is the first year, so prepare the luc data:

        ! open the luc file

        open(unit=15,file=lucdat(1:strlen(lucdat))//'.LUC')

        read (15,*) 
        read (15,*)
        read (15,*)

!       get first year of luc data
!       note we load the nfcancmx, not pfcancmx array. this is because this
!       nfcancmx value is passed to the pfcancmx array at the start of each simulation
!       year

        do i = 1, nltest
         if (.not. mosaic) then  !composite
           read (15,*) lucyr,(nfcancmxrow(i,1,j),j=1,icc)
         else                    !mosaic
           read (15,*) lucyr,(temparray(j),j=1,icc)
           do m = 1, nmtest-1 !as nmtest-1 = icc
             j = m
             nfcancmxrow(i,m,j) = temparray(m)
           enddo !m loop
         endif
        enddo !nltest

!       next update our luc data if either of the following conditions are met:
!       1) we are cycling the met data and the luc year we just read in is less
!       than the year we want to cycle over (assuming it is not defaulted to 
!       -9999) or,
!       2) we are not cycling over the met data so we just want to get the same
!       year of luc data as the met data we read in above in preparation for our
!       transient run.

        do while ((cyclemet .and. lucyr .lt. cylucyr             &
          .and. cylucyr .ne. -9999) .or. (.not. cyclemet .and.  &
          lucyr .lt. iyear))

!           get the luc data
            do i = 1, nltest
             if (.not. mosaic) then  !composite
               read (15,*,end=999) lucyr,(nfcancmxrow(i,1,j),j=1,icc)
             else                    !mosaic
               read (15,*,end=999) lucyr,(temparray(j),j=1,icc)
               do m = 1, nmtest-1 !nmtest-1 same as icc
                j = m
                nfcancmxrow(i,m,j) = temparray(m) 
               enddo !m loop
             endif
            enddo !nltest
        enddo  !while loop

!       get fcans for use by class using the nfcancmxs just read in
        k1=0
        do 997 j = 1, ican
          if(j.eq.1) then
            k1 = k1 + 1
          else
            k1 = k1 + nol2pfts(j-1)
          endif
          k2 = k1 + nol2pfts(j) - 1
          do 998 n = k1, k2
            do i = 1, nltest
             do m = 1, nmtest
              if (.not. mosaic) then !composite

               fcanrow(i,m,j)=fcanrow(i,m,j)+nfcancmxrow(i,m,n) 

              else if (mosaic .and. nfcancmxrow(i,m,n) .gt. seed) then
!              this tile has some plants so overwrite the seed fraction with
!              an actual fraction

!              note: the seed fraction has already been assigned in runclassctem
!              prior to entering this subroutine.
               farerow(i,m)=nfcancmxrow(i,m,n)    

              endif
             enddo
            enddo
998       continue
997     continue

!       (re)find the bare fraction for farerow(i,icc+1)
         if (mosaic) then
          do i = 1, nltest
          temp = 0.
           do m = 1, nmtest-1
            temp = temp + farerow(i,m)
           enddo

           farerow(i,nmtest) = 1.0 - temp

          enddo
         endif

          ! check that the bare fraction is possible (>0) and if not then
          ! reduce the other pfts proportionally to make a non-negative bare
          ! ground fraction.

          do i = 1, nltest
           temp = 0.0
           if (mosaic) then
            if (farerow(i,nmtest) < seed) then

             nfcancmxrow= adjust_luc_fracs(nlat,nmos,i,nmtest,nltest,icc, &
                                          nfcancmxrow,farerow(i,nmtest))

             do m = 1, nmtest
                n = m
                farerow(i,m)=nfcancmxrow(i,m,n)  
                temp = temp + farerow(i,m)  
             enddo

             farerow(i,nmtest) = 1.0 - temp

            endif !farerow<seed
           endif !mosaic
          enddo !nltest


!       assign the present pft fractions from those just read in
        do j = 1, icc
          do i = 1, nltest
           do m = 1, nmtest
            if (.not. mosaic) then  !composite 
             fcancmxrow(i,m,j)=nfcancmxrow(i,m,j)
             pfcancmxrow(i,m,j)=nfcancmxrow(i,m,j)
            else !mosaic
!            ensure that the fraction is >= seed
             pfcancmxrow(i,m,j)=max(seed,nfcancmxrow(i,m,j))
            endif
           enddo
          enddo
        enddo

!       back up one year in the luc file 

        do i = 1, nltest
           backspace(15)  
        enddo

return

999    continue
  
! end of the luc file is reached. close and tell main program to exit
        close(15)
        reach_eof = .true.

end subroutine initialize_luc

!=======================================================================

subroutine readin_luc(nlat,nmos,iyear,nmtest,nltest,mosaic,icc,lucyr, &
                      nfcancmxrow,reach_eof)

!           Canadian Terrestrial Ecosystem Model (CTEM) 
!                    Luc Annual Read-In Subroutine 
!
!     9  Jan. 2013  - this subroutine takes in luc information from
!     J. Melton       a luc file annually and adapts them for runclassctem
!		      

implicit none

! arguments

! inputs
integer, intent(in) :: nmos
integer, intent(in) :: nlat
integer, intent(in) :: iyear
logical, intent(in) :: mosaic
integer, intent(in) :: nmtest
integer, intent(in) :: nltest
integer, intent(in) :: icc

! updates
integer, intent(inout) :: lucyr
logical, intent(inout) :: reach_eof

! outputs
real, dimension(nlat,nmos,icc), intent(out) :: nfcancmxrow

! local variables
real, dimension(icc) :: temparray
real :: temp
integer :: j,m,i
integer :: k1,k2,n
real :: bare_ground_frac

!-------------------------

!it is subsequent years so read in and adjust the luc file info.

         do while (lucyr <= iyear) 
           do i = 1, nltest
            if (.not. mosaic) then  !composite
              read (15,*,end=999) lucyr,(nfcancmxrow(i,1,j),j=1,icc)
            else                    !mosaic
              read (15,*,end=999) lucyr,(temparray(j),j=1,icc)
              do m = 1, nmtest-1    !nmtest-1 same as icc
               j = m
               nfcancmxrow(i,m,j) = max(seed,temparray(m)) 
              enddo !m loop
            endif
           enddo !nltest
         enddo !lucyr<iyear
  
!       (re)find the bare fraction for farerow(i,icc+1)
         if (mosaic) then
          do i = 1, nltest
          temp = 0.0
           do m = 1, nmtest-1
            j = m
            temp = temp + nfcancmxrow(i,m,j)
           enddo

            bare_ground_frac = 1.0- temp

          enddo
          
          do i = 1, nltest
           if (bare_ground_frac < seed) then

             nfcancmxrow= adjust_luc_fracs(nlat,nmos,i,nmtest,nltest,icc, &
                                           nfcancmxrow,bare_ground_frac)

           endif !bare_ground_frac<seed
          enddo !nltest
         endif !mosaic

return

999 continue

! end of the luc file is reached. close and tell main program to exit
        close(15)
        reach_eof = .true.

end subroutine readin_luc

!=======================================================================

subroutine    luc(     icc,      ilg,      il1,       il2,  & 
                              ic, nol2pfts,    l2max,               & !1    
                        grclarea, pfcancmx, nfcancmx,      iday,    & !2    
                         todfrac,  yesfrac, interpol,               & !3    
!    ----------------------- inputs above this line -------------       
                         gleafmas, bleafmas, stemmass, rootmass,    & !4  
                         litrmass, soilcmas, vgbiomas, gavgltms,    & !5   
                         gavgscms,  fcancmx,   fcanmx,              & !6
!    ----------- updates above this line, outputs below ---------
                         lucemcom, lucltrin, lucsocin)                !7  
!
!     ----------------------------------------------------------------
!
!           Canadian Terrestrial Ecosystem Model (CTEM) 
!                       Land Use Change Subroutine 
!
!     18  Apr. 2013 - made it so that you will exit luc if the grid cell has
!     J. Melton       no actual luc in this timestep. removed some extraneous checks.
!
!     02  Jan. 2004 - this subroutine deals with the changes in the land
!     V. Arora        cover and estimates land use change (luc)
!                     related carbon emissions. based on whether there
!                     is conversion of forests/grassland to crop area,
!                     or croplands abandonment, this subroutine
!                     reallocates carbon from live vegetation to litter
!                     and soil c components. the decomposition from the
!                     litter and soil c pools thus implicitly models luc
!                     related carbon emissions. set of rules are
!                     followed to determine the fate of carbon that
!                     results from deforestation or replacement of
!                     grasslands by crops. 
!
!     ----------------------------------------------------------------
!     inputs
!
!     icc       - no of pfts for use by ctem, currently 9
!     ic        - no of pfts for use by class, currently 4
!     ilg      - no. of grid cells in latitude circle
!     il1, il2  - il1=1, il2=ilg
!     nol2pfts  - number of level 2 pfts
!     l2max     - maximum number of level 2 pfts
!     fcancmx   - max. fractional coverages of ctem's 9 pfts. 
!     grclarea  - gcm grid cell area, km2
!     pfcancmx  - previous max. fractional coverages of ctem's 9 pfts. 
!     nfcancmx  - next max. fractional coverages of ctem's 9 pfts. 
!     iday      - day of year
!     todfrac   - today's fractional coverage of all pfts 
!     yesfrac   - yesterday's fractional coverage of all pfts 
!     interpol  - if todfrac & yesfrac are provided then interpol must
!                 be set to false so that this subroutine doesn't do its
!                 own interpolation using pfcancmx and nfcancmx which
!                 are year end values     
!     updates

!     gleafmas  - green or live leaf mass in kg c/m2, for the 9 pfts
!     bleafmas  - brown or dead leaf mass in kg c/m2, for the 9 pfts
!     stemmass  - stem biomass in kg c/m2, for the 9 pfts
!     rootmass  - root biomass in kg c/m2, for the 9 pfts
!     litrmass  - litter mass in kg c/m2, for the 9 pfts + bare
!     soilcmas  - soil c mass in kg c/m2, for the 9 pfts + bare
!     vgbiomas  - grid averaged vegetation biomass, kg c/m2
!     gavgltms  - grid averaged litter mass, kg c/m2
!     gavgscms  - grid averaged soil c mass, kg c/m2

!     outputs

!     fcanmx   - fractional coverages of class 4 pfts (these are found based 
!                on new fcancmxs)
!     lucemcom - luc related carbon emission losses from combustion
!                u-mol co2/m2.sec
!     lucltrin - luc related input to litter pool, u-mol co2/m2.sec
!     lucsocin - luc related input to soil carbon pool, u-mol co2/m2.sec

!     ----------------------------------------------------------------    

      implicit none

      integer icc,  ilg, il1, il2, ic, i, j, k, m, n, kk, k1, k2, q 
      integer    iday, nol2pfts(ic),         l2max       
      integer fraciord(ilg,icc), treatind(ilg,icc),       bareiord(ilg)
      integer lrgstpft(1)

      logical  interpol, luctkplc(ilg) 
!      
      real  gleafmas(ilg,icc), bleafmas(ilg,icc),  stemmass(ilg,icc), & !1
            rootmass(ilg,icc),  fcancmx(ilg,icc),  pfcancmx(ilg,icc), & !2
                vgbiomas(ilg),             zero,  soilcmas(ilg,icc+1),& !3 
          litrmass(ilg,icc+1),     gavgltms(ilg),      gavgscms(ilg), & !4
            nfcancmx(ilg,icc),  fcancmy(ilg,icc),   todfrac(ilg,icc), & !5
             yesfrac(ilg,icc)

      real    fcanmx(ilg,ic),  delfrac(ilg,icc),          combust(3), & !1
                    paper(3),      furniture(3),   abvgmass(ilg,icc), & !2
                 bmasthrs(2),          grclarea,   combustc(ilg,icc),  & !3
             paperc(ilg,icc), furnturc(ilg,icc),  incrlitr(ilg,icc),  & !4
           incrsolc(ilg,icc),          tolrnce1,            tolrnce2, & !5
               chopedbm(ilg),           km2tom2

      real         redubmas1,              term,       barefrac(ilg), & 
                grsumcom(ilg),     grsumpap(ilg),      grsumfur(ilg), & !1
                grsumlit(ilg),     grsumsoc(ilg),      pbarefra(ilg), & !2
                grdencom(ilg),     grdenpap(ilg),      grdenfur(ilg), & !3
                grdenlit(ilg),     grdensoc(ilg),      totcmass(ilg), & !4
                totlmass(ilg),     totdmas1(ilg),      ntotcmas(ilg), & !5
                ntotlmas(ilg),     ntotdms1(ilg),      lucemcom(ilg), & !6
                pvgbioms(ilg),     pgavltms(ilg),      pgavscms(ilg), & !7
                    redubmas2,     lucltrin(ilg),       lucsocin(ilg),& !7
                totdmas2(ilg),     ntotdms2(ilg)


!     how much deforested/chopped off biomass is combusted
!     (these absolutely must add to 1.0!)
      data combust/0.15, 0.30, 0.45/ 

!     how much deforested/chopped off biomass goes into short term
!     storage such as paper
      data paper/0.70, 0.70, 0.55/

!     how much deforested/chopped off biomass goes into long term
!     storage such as furniture
      data furniture/0.15, 0.0, 0.0/

!     biomass thresholds for determining if deforested area is a forest,
!     a shrubland, or a bush kg c/m2
      data bmasthrs/4.0, 1.0/

      data zero/1.0e-20/

      data tolrnce1/0.50/  ! kg c, tolerance of total c balance 
      data tolrnce2/0.005/ ! kg c/m2, tolerance for total c density

      data km2tom2/1.0e+06/ ! changes from km2 to m2
!     ---------------------------------------------------------------

      if(icc.ne.9)                               call xit('luc',-1)  
      if(ic.ne.4)                                call xit('luc',-2)  

!     ------------------------------------------------------------------

!     find/use provided current and previous day's fractional coverage 

      if(interpol) then ! perform interpolation 
       do 110 j = 1, icc
        do 111 i = il1, il2

          delfrac(i,j)=nfcancmx(i,j)-pfcancmx(i,j) !change in fraction
          delfrac(i,j)=delfrac(i,j)/365.0
          fcancmx(i,j)=pfcancmx(i,j)+(real(iday)*delfrac(i,j)) !  current day
          fcancmy(i,j)=pfcancmx(i,j)+(real(iday-1)*delfrac(i,j)) ! previous day

          if( fcancmx(i,j).lt.0.0.and.abs(fcancmx(i,j)).lt.1.0e-05)then
            fcancmx(i,j)=0.0
          else if( fcancmx(i,j).lt.0.0.and. &
                abs(fcancmx(i,j)).ge.1.0e-05)then
            write(6,*)'fcancmx(',i,',',j,')=',fcancmx(i,j)
            write(6,*)'fractional coverage cannot be negative'
            call xit('luc',-4)
          endif

          if(fcancmy(i,j).lt.0.0.and.abs(fcancmy(i,j)).lt.1.0e-05)then    
            fcancmy(i,j)=0.0
          else if( fcancmy(i,j).lt.0.0.and. &
         abs(fcancmy(i,j)).ge.1.0e-05)then
            write(6,*)'fcancmy(',i,',',j,')=',fcancmy(i,j)
            write(6,*)'fractional coverage cannot be negative'
            call xit('luc',-5)
          endif

111     continue
110    continue
      else ! use provided values but still check
!            they are not -ve
       do 115 j = 1, icc
        do 116 i = il1, il2
          fcancmx(i,j) = todfrac(i,j)   
          fcancmy(i,j) = yesfrac(i,j)   

          if( fcancmx(i,j).lt.0.0.and.abs(fcancmx(i,j)).lt.1.0e-05)then
            fcancmx(i,j)=0.0
          else if( fcancmx(i,j).lt.0.0.and. &
         abs(fcancmx(i,j)).ge.1.0e-05)then
            write(6,*)'fcancmx(',i,',',j,')=',fcancmx(i,j)
            write(6,*)'fractional coverage cannot be negative'
            call xit('luc',-4)
          endif

          if(fcancmy(i,j).lt.0.0.and.abs(fcancmy(i,j)).lt.1.0e-05)then
            fcancmy(i,j)=0.0
          else if( fcancmy(i,j).lt.0.0.and. &
         abs(fcancmy(i,j)).ge.1.0e-05)then
            write(6,*)'fcancmy(',i,',',j,')=',fcancmy(i,j)
            write(6,*)'fractional coverage cannot be negative'
            call xit('luc',-5)
          endif

116     continue
115    continue
      endif

!     check if this year's fractional coverages have changed or not
!     for any pft

      do 150 i = il1, il2
        luctkplc(i)=.false.  ! did land use change take place for any pft
150   continue         ! in this grid cell

      do 200 j = 1, icc
        do 250 i = il1, il2
          if ( (abs(fcancmx(i,j)-fcancmy(i,j))).gt.zero ) then
            luctkplc(i)=.true. ! yes, luc did take place in this grid cell
          endif
250     continue
200   continue

!     only perform the rest of the subroutine if any luc is actually taking 
!     place, otherwise exit.
      do 255 q = il1, il2
        if (luctkplc(q)) then
!     -------------------------------------------------------------------
 
!     initialization

      do 260 j = 1, ic
        do 261 i = il1, il2
            fcanmx(i,j)=0.0 ! fractional coverage of class' pfts
261     continue      
260   continue      

      do 270 j = 1, icc
        do 271 i = il1, il2
          fraciord(i,j)=0   !fractional coverage increase or decrease
!                           !increase +1, decrease -1
          abvgmass(i,j)=0.0 !above-ground biomass
          treatind(i,j)=0   !treatment index for combust, paper, & furniture
          combustc(i,j)=0.0 !total carbon from deforestation- combustion
          paperc(i,j)=0.0   !total carbon from deforestation- paper
          furnturc(i,j)=0.0 !total carbon from deforestation- furniture

271     continue  
270   continue      

      do 280 i = il1, il2
        pvgbioms(i)=vgbiomas(i)  ! store grid average quantities in
        pgavltms(i)=gavgltms(i)  ! temporary arrays
        pgavscms(i)=gavgscms(i)

        vgbiomas(i)=0.0
        gavgltms(i)=0.0
        gavgscms(i)=0.0

        barefrac(i)=1.0          ! initialize bare fraction to 1.0
        pbarefra(i)=1.0          ! initialize previous years's bare fraction to 1.0

        grsumcom(i)=0.0          ! grid sum of combustion carbon for all
                                 ! pfts that are chopped
        grsumpap(i)=0.0          ! similarly for paper, 
        grsumfur(i)=0.0          ! furniture,
        grsumlit(i)=0.0          ! litter, and
        grsumsoc(i)=0.0          ! soil c

        grdencom(i)=0.0          ! grid averaged densities for combustion carbon,
        grdenpap(i)=0.0          ! paper,
        grdenfur(i)=0.0          ! furniture, 
        grdenlit(i)=0.0          ! litter, and
        grdensoc(i)=0.0          ! soil c

        totcmass(i)=0.0          ! total c mass (live+dead)
        totlmass(i)=0.0          ! total c mass (live)
        totdmas1(i)=0.0          ! total c mass (dead) litter
        totdmas2(i)=0.0          ! total c mass (dead) soil c

!       and the same 3 quantities as above after luc treatment

        ntotcmas(i)=0.0          ! total c mass (live+dead)
        ntotlmas(i)=0.0          ! total c mass (live)
        ntotdms1(i)=0.0          ! total c mass (dead) litter
        ntotdms2(i)=0.0          ! total c mass (dead) soil c

        lucemcom(i)=0.0          ! luc related combustion emission losses
        lucltrin(i)=0.0          ! luc related inputs to litter pool
        lucsocin(i)=0.0          ! luc related inputs to soil c pool

        bareiord(i)=0            ! bare fraction increases or decreases
        chopedbm(i)=0.0          ! chopped off biomass
280   continue

!     initialization ends

!     -------------------------------------------------------------------

!     if land use change has taken place then get fcanmxs for use by 
!     class based on the new fcancmxs

      k1=0
      do 300 j = 1, ic
        if(j.eq.1) then
          k1 = k1 + 1
        else
          k1 = k1 + nol2pfts(j-1)
        endif
        k2 = k1 + nol2pfts(j) - 1
        do 301 m = k1, k2
          do 302 i = il1, il2
            fcanmx(i,j)=fcanmx(i,j)+fcancmx(i,m)
            barefrac(i)=barefrac(i)-fcancmx(i,m)
302       continue
301     continue
300   continue

!     check if the interpol didn't mess up the barefrac. if so, take the
!     extra amount from the pft with the largest area. jm apr 24 2013.
      do 304 i = il1, il2
         if (barefrac(i).lt.0.0.and.abs(barefrac(i)).ge.1.0e-05) then
            lrgstpft = maxloc(fcancmx(i,1:icc))
            fcancmx(i,lrgstpft(1)) = fcancmx(i,lrgstpft(1)) + barefrac(i)
            barefrac(i) = 0.0
         endif
304   continue
!     find previous day's bare fraction using previous day's fcancmxs

      do 310 j = 1, icc
        do 311 i = il1, il2
          pbarefra(i)=pbarefra(i)-fcancmy(i,j)
311     continue
310   continue

!     check if the interpol didn't mess up the pbarefra. if so, take the
!     extra amount from the pft with the largest area. jm apr 24 2013.
      do 314 i = il1, il2
         if (pbarefra(i).lt.0.0.and.abs(pbarefra(i)).ge.1.0e-05) then
            lrgstpft = maxloc(fcancmy(i,1:icc))
            fcancmy(i,lrgstpft(1)) = fcancmy(i,lrgstpft(1)) + pbarefra(i)
            pbarefra(i) = 0.0
         endif
314   continue


!     based on sizes of 3 live pools and 2 dead pools we estimate the
!     total amount of c in each grid cell.

      do 320 j = 1, icc
        do 321 i = il1, il2
          totlmass(i)=totlmass(i)+ &
                     (fcancmy(i,j)*(gleafmas(i,j)+bleafmas(i,j)+ &
                      stemmass(i,j)+rootmass(i,j))*grclarea &
                      *km2tom2)
321     continue
320   continue

      do 340 j = 1, icc+1
        do 341 i = il1, il2
          if(j.lt.icc+1) then
            term = fcancmy(i,j)
          else if(j.eq.icc+1) then
            term = pbarefra(i)
          endif
          totdmas1(i)=totdmas1(i)+ &
                     (term*litrmass(i,j)*grclarea*km2tom2)
          totdmas2(i)=totdmas2(i)+ & 
                     (term*soilcmas(i,j)*grclarea*km2tom2)

341     continue
340   continue

      do 350 i = il1, il2
        totcmass(i)=totlmass(i)+totdmas1(i)+totdmas2(i)
350   continue

!     bare fractions cannot be negative

      do 440 i = il1, il2
        if( pbarefra(i).lt.0.0.and.abs(pbarefra(i)).lt.1.0e-05 )then
          pbarefra(i)=0.0
        else if(pbarefra(i).lt.0.0.and.abs(pbarefra(i)).ge.1.0e-05 )then
          write(6,*)'bare fractions cannot be negative'
          write(6,*)'prev. bare fraction(',i,')  =',pbarefra(i)
          call xit('luc',-7)
        endif

        if( barefrac(i).lt.0.0.and.abs(barefrac(i)).lt.1.0e-05 )then
          barefrac(i)=0.0
        else if(barefrac(i).lt.0.0.and.abs(barefrac(i)).ge.1.0e-05 )then
          write(6,*)'bare fractions cannot be negative'
          write(6,*)'bare fraction(',i,')  =',barefrac(i)
          call xit('luc',-8)
        endif
440   continue

!     find above ground biomass and treatment index for combust, paper,
!     and furniture

      k1=0
      do 500 j = 1, ic
        if(j.eq.1) then
          k1 = k1 + 1
        else
          k1 = k1 + nol2pfts(j-1)
        endif
        k2 = k1 + nol2pfts(j) - 1
        do 510 m = k1, k2
          do 520 i = il1, il2
            abvgmass(i,m)=gleafmas(i,m)+bleafmas(i,m)+stemmass(i,m)
            if(j.eq.1.or.j.eq.2) then  ! trees
              if(abvgmass(i,m).ge.bmasthrs(1)) then !forest
                treatind(i,m)=1
              else if (abvgmass(i,m).le.bmasthrs(2)) then !bush 
                treatind(i,m)=3
              else  !shrubland
                treatind(i,m)=2
              endif
            else                       !crops and grasses
              treatind(i,m)=3
            endif
520       continue
510     continue
500   continue

!     if treatment index is zero then something is wrong - we exit
!     jm edit -- isn't this impossible, so we don't need to check?
      do 530 j = 1, icc
        do 531 i = il1, il2
          if(treatind(i,j).eq.0)then
        write(6,*)'treatment index zero for grid cell ',i,' and pft',j
            call xit('luc',-9)
          endif
531     continue
530   continue

!     check if a pft's fractional cover is increasing or decreasing

      do 550 j = 1, icc 
        do 551 i = il1, il2
          if( ( fcancmx(i,j).gt.fcancmy(i,j)) .and. &
             (abs(fcancmy(i,j)-fcancmx(i,j)).gt.zero) ) then
              fraciord(i,j)=1  ! increasing
          else if( ( fcancmx(i,j).lt.fcancmy(i,j)) .and. &
                  (abs(fcancmy(i,j)-fcancmx(i,j)).gt.zero) ) then
              fraciord(i,j)=-1 ! decreasing
          endif
551     continue
550   continue

!     check if bare fraction increases of decreases

      do 560 i = il1, il2
        if( ( barefrac(i).gt.pbarefra(i)) .and. &
           (abs(pbarefra(i)-barefrac(i)).gt.zero) ) then
              bareiord(i)=1  !increasing
        else if ( ( barefrac(i).lt.pbarefra(i)) .and. &
                 (abs(pbarefra(i)-barefrac(i)).gt.zero) ) then
              bareiord(i)=-1 !decreasing
        endif
560   continue

      
!     if the fractional coverage of pfts increases then spread their
!     live & dead biomass uniformly over the new fraction. this 
!     effectively reduces their per m2 c density. 
 
      do 570 j = 1, icc
        do 571 i = il1, il2
          if(fraciord(i,j).eq.1)then
            term = fcancmy(i,j)/fcancmx(i,j)
            gleafmas(i,j)=gleafmas(i,j)*term
            bleafmas(i,j)=bleafmas(i,j)*term
            stemmass(i,j)=stemmass(i,j)*term
            rootmass(i,j)=rootmass(i,j)*term
            litrmass(i,j)=litrmass(i,j)*term
            soilcmas(i,j)=soilcmas(i,j)*term
          endif 
571     continue
570   continue

!     if bare fraction increases then spread its litter and soil c
!     uniformly over the increased fraction

      do 580 i = il1, il2
        if(bareiord(i).eq.1)then
          term = pbarefra(i)/barefrac(i)
          litrmass(i,icc+1)=litrmass(i,icc+1)*term
          soilcmas(i,icc+1)=soilcmas(i,icc+1)*term
        endif
580   continue

!     if any of the pfts fractional coverage decreases, then we chop the
!     aboveground biomass and treat it according to our rules (burn it,
!     and convert it into paper and furniture). the below ground live
!     biomass and litter of this pfts gets assimilated into litter of
!     all pfts (uniformly spread over the whole grid cell), and soil c
!     from the chopped off fraction of this pft, gets assimilated into
!     soil c of all existing pfts as well.

      k1=0
      do 600 j = 1, ic
        if(j.eq.1) then
          k1 = k1 + 1
        else
          k1 = k1 + nol2pfts(j-1)
        endif
        k2 = k1 + nol2pfts(j) - 1
        do 610 m = k1, k2
          do 620 i = il1, il2

            if(fraciord(i,m).eq.-1)then

!             chop off above ground biomass 
              redubmas1=(fcancmy(i,m)-fcancmx(i,m))*grclarea &
                       *abvgmass(i,m)*km2tom2

              if(redubmas1.lt.0.0)then
                write(6,*)'redubmas1 less than zero'
                write(6,*)'fcancmy-fcancmx = ',  &
                          fcancmy(i,m)-fcancmx(i,m)
                write(6,*)'grid cell = ',i,' pft = ',m
                call xit('luc',-10)
              endif

!             rootmass needs to be chopped as well and all of it goes to
!             the litter/paper pool

              redubmas2=(fcancmy(i,m)-fcancmx(i,m))*grclarea  &
                       *rootmass(i,m)*km2tom2

!             keep adding chopped off biomass for each pft to get the total
!             for a grid cell for diagnostics

              chopedbm(i)=chopedbm(i) + redubmas1 + redubmas2 

!             find what's burnt, and what's converted to paper &
!             furniture
              combustc(i,m)=combust(treatind(i,m))*redubmas1
              paperc(i,m)=paper(treatind(i,m))*redubmas1 + redubmas2
              furnturc(i,m)=furniture(treatind(i,m))*redubmas1

!             keep adding all this for a given grid cell
              grsumcom(i)=grsumcom(i)+combustc(i,m)
              grsumpap(i)=grsumpap(i)+paperc(i,m)
              grsumfur(i)=grsumfur(i)+furnturc(i,m)

!             litter from the chopped off fraction of the chopped
!             off pft needs to be assimilated, and so does soil c from
!             the chopped off fraction of the chopped pft

              redubmas1=(fcancmy(i,m)-fcancmx(i,m))*grclarea &
                       *litrmass(i,m)*km2tom2
              incrlitr(i,m)=redubmas1

              redubmas1=(fcancmy(i,m)-fcancmx(i,m))*grclarea &
                       *soilcmas(i,m)*km2tom2
              incrsolc(i,m)=redubmas1

              grsumlit(i)=grsumlit(i)+incrlitr(i,m)
              grsumsoc(i)=grsumsoc(i)+incrsolc(i,m)
            endif

620       continue
610     continue
600   continue

!     if bare fraction decreases then chop off the litter and soil c
!     from the decreased fraction and add it to grsumlit & grsumsoc
!     for spreading over the whole grid cell

      do 630 i = il1, il2
        if(bareiord(i).eq.-1)then

          redubmas1=(pbarefra(i)-barefrac(i))*grclarea &
                   *litrmass(i,icc+1)*km2tom2

          redubmas2=(pbarefra(i)-barefrac(i))*grclarea &
                   *soilcmas(i,icc+1)*km2tom2

          grsumlit(i)=grsumlit(i)+redubmas1
          grsumsoc(i)=grsumsoc(i)+redubmas2
        endif
630   continue

!     calculate if the chopped off biomass equals the sum of grsumcom(i),
!     grsumpap(i) & grsumfur(i)

      do 632 i = il1, il2
       if( abs(chopedbm(i)-grsumcom(i)-grsumpap(i)-grsumfur(i)).gt. &
          tolrnce1 ) then
           write(6,*)'at grid cell = ',i
           write(6,*)'chopped biomass does not equals sum of total'   
           write(6,*)'luc related emissions'
           write(6,*)'chopedbm(i) = ',chopedbm(i)
           write(6,*)'grsumcom(i) = ',grsumcom(i)
           write(6,*)'grsumpap(i) = ',grsumpap(i)
           write(6,*)'grsumfur(i) = ',grsumfur(i)
           write(6,*)'sum of grsumcom, grsumpap, grsumfur(i) = ', &
            grsumcom(i)+grsumpap(i)+grsumfur(i) 
           call xit('luc',-11)
       endif
632   continue     

!     spread chopped off stuff uniformly over the litter and soil c
!     pools of all existing pfts, including the bare fraction.

!     convert the available c into density 

      do 640 i = il1, il2
        grdencom(i)=grsumcom(i)/(grclarea*km2tom2)
        grdenpap(i)=grsumpap(i)/(grclarea*km2tom2)
        grdenfur(i)=grsumfur(i)/(grclarea*km2tom2)
        grdenlit(i)=grsumlit(i)/(grclarea*km2tom2)
        grdensoc(i)=grsumsoc(i)/(grclarea*km2tom2)

640   continue

      do 650 j = 1, icc
        do 651 i = il1, il2
          if(fcancmx(i,j).gt.zero)then
            litrmass(i,j)=litrmass(i,j)+grdenpap(i)+grdenlit(i)
            soilcmas(i,j)=soilcmas(i,j)+grdenfur(i)+grdensoc(i)
          else
            gleafmas(i,j)=0.0
            bleafmas(i,j)=0.0
            stemmass(i,j)=0.0
            rootmass(i,j)=0.0
            litrmass(i,j)=0.0
            soilcmas(i,j)=0.0
          endif
651     continue 
650   continue
 
      do 660 i = il1, il2
        if(barefrac(i).gt.zero)then
          litrmass(i,icc+1)=litrmass(i,icc+1)+grdenpap(i)+grdenlit(i)
          soilcmas(i,icc+1)=soilcmas(i,icc+1)+grdenfur(i)+grdensoc(i)
        else
          litrmass(i,icc+1)=0.0
          soilcmas(i,icc+1)=0.0
        endif
660   continue

!     the combusted c is used to find the c flux that we can release
!     into the atmosphere.

      do 670 i = il1, il2
        lucemcom(i)=grdencom(i)     ! this is flux in kg c/m2.day that
!                                   ! will be emitted 
!       lucltrin(i)=grdenpap(i)+grdenlit(i) ! flux in kg c/m2.day
!       lucsocin(i)=grdenfur(i)+grdensoc(i) ! flux in kg c/m2.day

        lucltrin(i)=grdenpap(i) ! flux in kg c/m2.day
        lucsocin(i)=grdenfur(i) ! flux in kg c/m2.day

!       convert all land use change fluxes to u-mol co2-c/m2.sec
        lucemcom(i)=lucemcom(i)*963.62
        lucltrin(i)=lucltrin(i)*963.62
        lucsocin(i)=lucsocin(i)*963.62

670   continue

!     and finally we see if the total amount of carbon is conserved

      do 700 j = 1, icc
        do 701 i = il1, il2
          ntotlmas(i)=ntotlmas(i)+ &
                     (fcancmx(i,j)*(gleafmas(i,j)+bleafmas(i,j)+&
                     stemmass(i,j)+rootmass(i,j))*grclarea &
                     *km2tom2)
701     continue
700   continue

      do 710 j = 1, icc+1
        do 711 i = il1, il2
          if(j.lt.icc+1) then
            term = fcancmx(i,j)
          else if(j.eq.icc+1) then
            term = barefrac(i)
          endif
          ntotdms1(i)=ntotdms1(i)+ &
                     (term*litrmass(i,j)*grclarea*km2tom2) 
          ntotdms2(i)=ntotdms2(i)+ & 
                     (term*soilcmas(i,j)*grclarea*km2tom2)
711     continue
710   continue

      do 720 i = il1, il2
        ntotcmas(i)=ntotlmas(i)+ntotdms1(i)+ntotdms2(i)
720   continue

!     total litter mass (before + input from chopped off biomass)
!     and after must be same

      do 722 i = il1, il2
        if( abs(totdmas1(i)+grsumpap(i)-ntotdms1(i)).gt.tolrnce1 )then   
           write(6,*)'at grid cell = ',i
           write(6,*)'total litter carbon does not balance after luc'
           write(6,*)'totdmas1(i) = ',totdmas1(i)
           write(6,*)'grsumpap(i) = ',grsumpap(i)
           write(6,*)'totdmas1(i) + grsumpap(i) = ', &
                     totdmas1(i) + grsumpap(i)
           write(6,*)'ntotdms1(i) = ',ntotdms1(i)
           call xit('luc',-12)
        endif
722   continue

!     for conservation totcmass(i) must be equal to ntotcmas(i) plus
!     combustion carbon losses

      do 730 i = il1, il2
        if( abs(totcmass(i)-ntotcmas(i)-grsumcom(i)).gt.tolrnce1)then
           write(6,*)'at grid cell = ',i
           write(6,*)'total carbon does not balance after luc'
           write(6,*)'totcmass(i) = ',totcmass(i)
           write(6,*)'ntotcmas(i) = ',ntotcmas(i)
           write(6,*)'grsumcom(i) = ',grsumcom(i)
           call xit('luc',-13)
        endif
730   continue

!     update grid averaged vegetation biomass, and litter and soil c
!     densities

      do 750 j = 1, icc
        do 751 i = il1, il2
          vgbiomas(i)=vgbiomas(i)+fcancmx(i,j)*(gleafmas(i,j)+ &
                     bleafmas(i,j)+stemmass(i,j)+rootmass(i,j))
          gavgltms(i)=gavgltms(i)+fcancmx(i,j)*litrmass(i,j)
          gavgscms(i)=gavgscms(i)+fcancmx(i,j)*soilcmas(i,j)
751     continue
750   continue

      do 760 i = il1, il2
        gavgltms(i)=gavgltms(i)+( barefrac(i)*litrmass(i,icc+1) )
        gavgscms(i)=gavgscms(i)+( barefrac(i)*soilcmas(i,icc+1) )
760   continue

!     just like total amount of carbon must balance, the grid averagred
!     densities must also balance

      do 780 i = il1, il2
       if( abs(pvgbioms(i)+pgavltms(i)+pgavscms(i)- &
              vgbiomas(i)-gavgltms(i)-gavgscms(i)- &
              grdencom(i)).gt.tolrnce2 ) then
           write(6,*)'iday = ',iday
           write(6,*)'at grid cell = ',i
           write(6,*)'pbarefra(i) = ',pbarefra(i)
           write(6,*)'barefrac(i) = ',barefrac(i)
           write(6,*)'pfcancmx(i,j) = ',(pfcancmx(i,j),j=1,icc)
           write(6,*)'nfcancmx(i,j) = ',(nfcancmx(i,j),j=1,icc)
           write(6,*)'total carbon density does not balance after luc'
           write(6,*)'pvgbioms(i) = ',pvgbioms(i)
           write(6,*)'pgavltms(i) = ',pgavltms(i)
           write(6,*)'pgavscms(i) = ',pgavscms(i)
           write(6,*)'vgbiomas(i) = ',vgbiomas(i)
           write(6,*)'gavgltms(i) = ',gavgltms(i)
           write(6,*)'gavgscms(i) = ',gavgscms(i)
           write(6,*)'grdencom(i) = ',grdencom(i)
           write(6,*)'pvgbioms + pgavltms + pgavscms = ', &
                   (pvgbioms(i)+pgavltms(i)+pgavscms(i))
           write(6,*)'vgbiomas + gavgltms + gavgscms + grdencom = ', &
         (vgbiomas(i)+gavgltms(i)+gavgscms(i)+ grdencom(i))
         write(6,*)'diff = ',abs((pvgbioms(i)+pgavltms(i)+pgavscms(i)) &
          -(vgbiomas(i)+gavgltms(i)+gavgscms(i)+ grdencom(i)))
           write(6,*)'tolrnce2 = ',tolrnce2
           call xit('luc',-14)
       endif
780   continue

      do 800 j = 1, icc
        do 810 i = il1, il2
          if(iday.eq.365)then
            pfcancmx(i,j)=nfcancmx(i,j)
          endif
810     continue
800   continue

       endif  ! loop to check if any luc took place. 
255   continue

      return 

end subroutine luc

!=======================================================================

function adjust_luc_fracs(nlat,nmos,i,nmtest,nltest,icc,nfcancmxrow, &
                          bare_ground_frac) result(outnfcrow)

! this function adjusts the amount of each pft to ensure that the fraction
! of gridcell bare ground is >0.

! j. melton, jan 11 2013

  implicit none

! arguments:
integer, intent(in) :: nmos
integer, intent(in) :: nlat
integer, intent(in) :: i
integer, intent(in) :: nmtest
integer, intent(in) :: nltest
integer, intent(in) :: icc
real, dimension(nlat,nmos,icc), intent(in) :: nfcancmxrow
real, intent(in) :: bare_ground_frac

real, dimension(nlat,nmos,icc) :: outnfcrow

! local variables:
integer :: m, j
real, dimension(icc) :: frac_abv_seed
real :: needed_bare,reduce_per_pft,tot

!-------------------------
tot = 0.0

! find the amount of needed space in the other pfts
! need a minimum bare area of seed.
needed_bare=(-bare_ground_frac) + seed

! get the proportionate amounts of each pft above the seed lower limit
do j = 1,nmtest-1
  m = j
  frac_abv_seed(j)=max(0.0,nfcancmxrow(i,m,j)-seed)
  tot = tot + frac_abv_seed(j)
enddo

! add in the bare ground min fraction of seed.
  tot = tot + seed

! now reduce the pfts proportional to their fractional area to make
! the bare ground fraction be the seed amount and no other pft be less than
! seed
do j = 1,nmtest-1
  m = j
  outnfcrow(i,m,j)=max(seed,nfcancmxrow(i,m,j) - needed_bare * &
                       frac_abv_seed(j) / tot)
enddo

end function adjust_luc_fracs

end module
