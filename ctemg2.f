      subroutine ctemg2(fcancmxgat,     rmatcgat,      zolncgat,
     1      paicgat,
     1      ailcgat,     ailcggat,    cmasvegcgat,     slaicgat,
     2      ailcgsgat,   fcancsgat,   fcancgat,     rmatctemgat,
     3      co2concgat,  co2i1cggat,  co2i1csgat,   co2i2cggat, 
     4      co2i2csgat,  xdiffusgat,  slaigat,      cfluxcggat, 
     5      cfluxcsgat,  ancsveggat,  ancgveggat,   rmlcsveggat,
     6      rmlcgveggat, canresgat,   sdepgat,
     7      sandgat,     claygat,     orgmgat,
     8      anveggat,    rmlveggat,   tcanoaccgat_m,tbaraccgat_m,
     9      uvaccgat_m,  vvaccgat_m,  mlightnggat,  prbfrhucgat,
     a      extnprobgat, stdalngat,   pfcancmxgat,  nfcancmxgat,
     b      stemmassgat, rootmassgat, litrmassgat,  gleafmasgat,
     c      bleafmasgat, soilcmasgat, ailcbgat,     flhrlossgat,
     d      pandaysgat,  lfstatusgat, grwtheffgat,  lystmmasgat,
     e      lyrotmasgat, tymaxlaigat, vgbiomasgat,  gavgltmsgat,
     f      stmhrlosgat, bmasveggat,  colddaysgat,  rothrlosgat,
     g      alvsctmgat,  alirctmgat,  gavglaigat,   nppgat,
     h      nepgat,      hetroresgat, autoresgat,   soilrespgat,
     i      rmgat,       rggat,       nbpgat,       litresgat,
     j      socresgat,   gppgat,      dstcemlsgat,  litrfallgat,
     k      humiftrsgat, veghghtgat,  rootdpthgat,  rmlgat,
     l      rmsgat,      rmrgat,      tltrleafgat,  tltrstemgat,
     m      tltrrootgat, leaflitrgat, roottempgat,  afrleafgat,
     n      afrstemgat,  afrrootgat,  wtstatusgat,  ltstatusgat,
     o      burnfracgat, probfiregat, lucemcomgat,  lucltringat,
     p      lucsocingat, nppveggat,   dstcemls3gat,
     q      faregat,     gavgscmsgat, rmlvegaccgat, pftexistgat,
     r      rmsveggat,   rmrveggat,   rgveggat,    vgbiomas_veggat,
     s      gppveggat,   nepveggat,   ailcmingat,   ailcmaxgat,
     t      emit_co2gat,  emit_cogat, emit_ch4gat,  emit_nmhcgat,
     u      emit_h2gat,   emit_noxgat,emit_n2ogat,  emit_pm25gat,
     v      emit_tpmgat,  emit_tcgat, emit_ocgat,   emit_bcgat,
     &      nbpveggat, hetroresveggat, autoresveggat,litresveggat,
     &      soilcresveggat,burnveggat,
c
     r      ilmos,       jlmos,       iwmos,        jwmos,
     s      nml,    fcancmxrow,  rmatcrow,    zolncrow,     paicrow,
     v      ailcrow,     ailcgrow,    cmasvegcrow,  slaicrow,
     w      ailcgsrow,   fcancsrow,   fcancrow,     rmatctemrow,
     x      co2concrow,  co2i1cgrow,  co2i1csrow,   co2i2cgrow,
     y      co2i2csrow,  xdiffus,     slairow,      cfluxcgrow,
     z      cfluxcsrow,  ancsvegrow,  ancgvegrow,   rmlcsvegrow,
     1      rmlcgvegrow, canresrow,   sdeprow,
     2      sandrow,     clayrow,     orgmrow,
     3      anvegrow,    rmlvegrow,   tcanoaccrow_m,tbaraccrow_m,
     4      uvaccrow_m,  vvaccrow_m,  mlightnggrd,  prbfrhucgrd,
     5      extnprobgrd, stdalngrd,   pfcancmxrow,  nfcancmxrow,
     6      stemmassrow, rootmassrow, litrmassrow,  gleafmasrow,
     7      bleafmasrow, soilcmasrow, ailcbrow,     flhrlossrow,
     8      pandaysrow,  lfstatusrow, grwtheffrow,  lystmmasrow,
     9      lyrotmasrow, tymaxlairow, vgbiomasrow,  gavgltmsrow,
     a      stmhrlosrow, bmasvegrow,  colddaysrow,  rothrlosrow,
     b      alvsctmrow,  alirctmrow,  gavglairow,   npprow,
     c      neprow,      hetroresrow, autoresrow,   soilresprow,
     d      rmrow,       rgrow,       nbprow,       litresrow,
     e      socresrow,   gpprow,      dstcemlsrow,  litrfallrow,
     f      humiftrsrow, veghghtrow,  rootdpthrow,  rmlrow,
     g      rmsrow,      rmrrow,      tltrleafrow,  tltrstemrow,
     h      tltrrootrow, leaflitrrow, roottemprow,  afrleafrow,
     i      afrstemrow,  afrrootrow,  wtstatusrow,  ltstatusrow,
     j      burnfracrow, probfirerow, lucemcomrow,  lucltrinrow,
     k      lucsocinrow, nppvegrow,   dstcemls3row,
     l      farerow,     gavgscmsrow, rmlvegaccrow, pftexistrow,
     &      rmsvegrow,   rmrvegrow,   rgvegrow,   vgbiomas_vegrow,
     &      gppvegrow,   nepvegrow,  ailcminrow,   ailcmaxrow,
     &      emit_co2row,  emit_corow, emit_ch4row,  emit_nmhcrow,
     &      emit_h2row,   emit_noxrow,emit_n2orow,  emit_pm25row,
     &      emit_tpmrow,  emit_tcrow, emit_ocrow,   emit_bcrow,
     &      nbpvegrow,   hetroresvegrow, autoresvegrow,litresvegrow,
     &      soilcresvegrow, burnvegrow )
c
c              Canadian Terrestrial Ecosystem Model (CTEM)
c
C     July 12 2013    Bring in the ctem params use statement
c     J. Melton
 
c     July 28 2009    Gather operation on CTEM variables.
c     Rong Li
c 
      use ctem_params,        only : nlat, nmos, ilg, ignd, ican, icp1,
     1                               icc,iccp1

      implicit none
c
c     * integer constants.
c
      integer  nml, k, l, m
c
c     * gather-scatter index arrays.
c
      integer  ilmos (ilg),  jlmos  (ilg),  iwmos  (ilg),  jwmos (ilg)
c
c
      real  fcancmxgat(ilg,icc),      rmatcgat(ilg,ican,ignd),
     1      zolncgat(ilg,ican),       paicgat(ilg,ican),
     2      ailcgat(ilg,ican),        ailcggat(ilg,icc),
     3      cmasvegcgat(ilg,ican),    slaicgat(ilg,ican),
     4      ailcgsgat(ilg,icc),       fcancsgat(ilg,icc),
     5      fcancgat(ilg,icc),        rmatctemgat(ilg,icc,ignd),
     6      co2concgat(ilg),          co2i1cggat(ilg,icc),
     7      co2i1csgat(ilg,icc),      co2i2cggat(ilg,icc),
     8      co2i2csgat(ilg,icc),      xdiffusgat(ilg),
     9      slaigat(ilg,icc),         cfluxcggat(ilg),
     a      cfluxcsgat(ilg),          ancsveggat(ilg,icc),
     b      ancgveggat(ilg,icc),      rmlcsveggat(ilg,icc),
     c      rmlcgveggat(ilg,icc),     canresgat(ilg),
     d      sdepgat(ilg)
c
      real    sandgat(ilg,ignd),      claygat(ilg,ignd), 
     1        orgmgat(ilg,ignd)
c
      real    anveggat(ilg,icc),      rmlveggat(ilg,icc)
c
      real   tcanoaccgat_m(ilg),
     1          uvaccgat_m(ilg),      vvaccgat_m(ilg)
c
      real   mlightnggat(ilg,12)  !12month
c
      real  prbfrhucgat(ilg),         extnprobgat(ilg),
     1      tbaraccgat_m(ilg,ignd),
     2      pfcancmxgat(ilg,icc),     nfcancmxgat(ilg,icc),
     3      stemmassgat(ilg,icc),     rootmassgat(ilg,icc),  
     4      litrmassgat(ilg,icc+1),   gleafmasgat(ilg,icc),
     5      bleafmasgat(ilg,icc),     soilcmasgat(ilg,icc+1),
     6      ailcbgat(ilg,icc),        flhrlossgat(ilg,icc)
c
      integer  pandaysgat(ilg,icc),   lfstatusgat(ilg,icc),
     1         stdalngat(ilg),       colddaysgat(ilg,2)
     
      logical pftexistgat(ilg,icc)
c
      real  grwtheffgat(ilg,icc),    lystmmasgat(ilg,icc),
     9      lyrotmasgat(ilg,icc),    tymaxlaigat(ilg,icc),
     a      vgbiomasgat(ilg),        gavgltmsgat(ilg),
     b      stmhrlosgat(ilg,icc),    bmasveggat(ilg,icc),
     c      rothrlosgat(ilg,icc),
     d      alvsctmgat(ilg,ican),    alirctmgat(ilg,ican),
     e      gavglaigat(ilg)
c
      real  nppgat(ilg),            nepgat(ilg),
     1      hetroresgat(ilg),       autoresgat(ilg),
     2      soilrespgat(ilg),       rmgat(ilg),     rggat(ilg),
     3      nbpgat(ilg),            litresgat(ilg),
     4      socresgat(ilg),         gppgat(ilg),
     5      dstcemlsgat(ilg),       litrfallgat(ilg),
     6      humiftrsgat(ilg),       veghghtgat(ilg,icc),
     7      rootdpthgat(ilg,icc),   rmlgat(ilg),
     8      rmsgat(ilg),            rmrgat(ilg),
     9      tltrleafgat(ilg,icc),   tltrstemgat(ilg,icc),
     a      tltrrootgat(ilg,icc),   leaflitrgat(ilg,icc),
     b      roottempgat(ilg,icc),   afrleafgat(ilg,icc),
     c      afrstemgat(ilg,icc),    afrrootgat(ilg,icc),
     d      wtstatusgat(ilg,icc),   ltstatusgat(ilg,icc),
     e      burnfracgat(ilg),       probfiregat(ilg),
     f      lucemcomgat(ilg),       lucltringat(ilg),
     g      lucsocingat(ilg),       nppveggat(ilg,icc),
     h      dstcemls3gat(ilg)
c
c      fire emission variables
       real emit_co2gat(ilg,icc),  emit_cogat(ilg,icc), 
     1      emit_ch4gat(ilg,icc),  emit_nmhcgat(ilg,icc),
     2      emit_h2gat(ilg,icc),   emit_noxgat(ilg,icc),
     3      emit_n2ogat(ilg,icc),  emit_pm25gat(ilg,icc),
     4      emit_tpmgat(ilg,icc),  emit_tcgat(ilg,icc),
     5      emit_ocgat(ilg,icc),   emit_bcgat(ilg,icc),
     6      burnveggat(ilg,icc)
c
      real  faregat(ilg)
      real  gavgscmsgat(ilg)
c
      real rmlvegaccgat(ilg,icc),       rmsveggat(ilg,icc),
     1      rmrveggat(ilg,icc),         rgveggat(ilg,icc),
     2      ailcmingat(ilg,icc),        ailcmaxgat(ilg,icc)
c
      real vgbiomas_veggat(ilg,icc)
c
      real gppveggat(ilg,icc),        nepveggat(ilg,icc),
     1     nbpveggat(ilg,icc),    hetroresveggat(ilg,iccp1),
     2      autoresveggat(ilg,icc),litresveggat(ilg,iccp1),
     3      soilcresveggat(ilg,iccp1) 
c
      real  fcancmxrow(nlat,nmos,icc),   rmatcrow(nlat,nmos,ican,ignd),
     1      zolncrow(nlat,nmos,ican),    paicrow(nlat,nmos,ican),
     2      ailcrow(nlat,nmos,ican),     ailcgrow(nlat,nmos,icc),
     3      cmasvegcrow(nlat,nmos,ican), slaicrow(nlat,nmos,ican),
     4      ailcgsrow(nlat,nmos,icc),    fcancsrow(nlat,nmos,icc),
     5      fcancrow(nlat,nmos,icc),    rmatctemrow(nlat,nmos,icc,ignd),
     6      co2concrow(nlat,nmos),       co2i1cgrow(nlat,nmos,icc),
     7      co2i1csrow(nlat,nmos,icc),   co2i2cgrow(nlat,nmos,icc),
     8      co2i2csrow(nlat,nmos,icc),   xdiffus(nlat),
     9      slairow(nlat,nmos,icc),      cfluxcgrow(nlat,nmos),
     a      cfluxcsrow(nlat,nmos),       ancsvegrow(nlat,nmos,icc),
     b      ancgvegrow(nlat,nmos,icc),   rmlcsvegrow(nlat,nmos,icc),
     c      rmlcgvegrow(nlat,nmos,icc),  canresrow(nlat,nmos),
     d      sdeprow(nlat,nmos)
c
      real    sandrow(nlat,nmos,ignd),   clayrow(nlat,nmos,ignd), 
     1        orgmrow(nlat,nmos,ignd)
c
      real  anvegrow(nlat,nmos,icc),     rmlvegrow(nlat,nmos,icc)
c
      real  tcanoaccrow_m(nlat,nmos),
     1      uvaccrow_m(nlat,nmos),       vvaccrow_m(nlat,nmos)
c
      real  mlightnggrd(nlat,12)  !12month
c
      real  prbfrhucgrd(nlat),            extnprobgrd(nlat),
     1      tbaraccrow_m(nlat,nmos,ignd),
     2      pfcancmxrow(nlat,nmos,icc),   nfcancmxrow(nlat,nmos,icc),
     3      stemmassrow(nlat,nmos,icc),   rootmassrow(nlat,nmos,icc),
     4      litrmassrow(nlat,nmos,icc+1), gleafmasrow(nlat,nmos,icc),
     5      bleafmasrow(nlat,nmos,icc),   soilcmasrow(nlat,nmos,icc+1),
     6      ailcbrow(nlat,nmos,icc),      flhrlossrow(nlat,nmos,icc)
c
      integer pandaysrow(nlat,nmos,icc),  lfstatusrow(nlat,nmos,icc),
     1      stdalngrd(nlat),               colddaysrow(nlat,nmos,2)

      logical pftexistrow(nlat,nmos,icc)
c
      real  grwtheffrow(nlat,nmos,icc),    lystmmasrow(nlat,nmos,icc),
     9      lyrotmasrow(nlat,nmos,icc),    tymaxlairow(nlat,nmos,icc),
     a      vgbiomasrow(nlat,nmos),        gavgltmsrow(nlat,nmos),
     b      stmhrlosrow(nlat,nmos,icc),    bmasvegrow(nlat,nmos,icc),
     c      rothrlosrow(nlat,nmos,icc),
     d      alvsctmrow(nlat,nmos,ican),    alirctmrow(nlat,nmos,ican),
     e      gavglairow(nlat,nmos)
c
      real  npprow(nlat,nmos),             neprow(nlat,nmos),
     1      hetroresrow(nlat,nmos),        autoresrow(nlat,nmos),
     2      soilresprow(nlat,nmos),        rmrow(nlat,nmos),
     2      rgrow(nlat,nmos),
     3      nbprow(nlat,nmos),             litresrow(nlat,nmos),
     4      socresrow(nlat,nmos),          gpprow(nlat,nmos),
     5      dstcemlsrow(nlat,nmos),        litrfallrow(nlat,nmos),
     6      humiftrsrow(nlat,nmos),        veghghtrow(nlat,nmos,icc),
     7      rootdpthrow(nlat,nmos,icc),    rmlrow(nlat,nmos),
     8      rmsrow(nlat,nmos),             rmrrow(nlat,nmos),
     9      tltrleafrow(nlat,nmos,icc),    tltrstemrow(nlat,nmos,icc),
     a      tltrrootrow(nlat,nmos,icc),    leaflitrrow(nlat,nmos,icc),
     b      roottemprow(nlat,nmos,icc),    afrleafrow(nlat,nmos,icc),
     c      afrstemrow(nlat,nmos,icc),     afrrootrow(nlat,nmos,icc),
     d      wtstatusrow(nlat,nmos,icc),    ltstatusrow(nlat,nmos,icc),
     e      burnfracrow(nlat,nmos),        probfirerow(nlat,nmos),
     f      lucemcomrow(nlat,nmos),        lucltrinrow(nlat,nmos),
     g      lucsocinrow(nlat,nmos),        nppvegrow(nlat,nmos,icc),
     h      dstcemls3row(nlat,nmos)
c
c     fire variables
      real emit_co2row(nlat,nmos,icc),  emit_corow(nlat,nmos,icc), 
     1     emit_ch4row(nlat,nmos,icc),  emit_nmhcrow(nlat,nmos,icc),
     2     emit_h2row(nlat,nmos,icc),   emit_noxrow(nlat,nmos,icc),
     3     emit_n2orow(nlat,nmos,icc),  emit_pm25row(nlat,nmos,icc),
     4     emit_tpmrow(nlat,nmos,icc),  emit_tcrow(nlat,nmos,icc),
     5     emit_ocrow(nlat,nmos,icc),   emit_bcrow(nlat,nmos,icc),
     6     burnvegrow(nlat,nmos,icc)
c
      real  farerow(nlat,nmos)
c
      real  gavgscmsrow(nlat,nmos)
c
      real rmlvegaccrow(nlat,nmos,icc), rmsvegrow(nlat,nmos,icc),
     1      rmrvegrow(nlat,nmos,icc),   rgvegrow(nlat,nmos,icc),
     2      ailcminrow(nlat,nmos,icc),  ailcmaxrow(nlat,nmos,icc)
c
      real vgbiomas_vegrow(nlat,nmos,icc)
c
      real gppvegrow(nlat,nmos,icc),    nepvegrow(nlat,nmos,icc),
     1     nbpvegrow(nlat,nmos,icc), hetroresvegrow(nlat,nmos,iccp1),
     2      autoresvegrow(nlat,nmos,icc),litresvegrow(nlat,nmos,iccp1),
     3      soilcresvegrow(nlat,nmos,iccp1)
c
c----------------------------------------------------------------------
      do 100 k=1,nml
          sdepgat(k)      = sdeprow(ilmos(k),jlmos(k))
          co2concgat(k)   = co2concrow(ilmos(k),jlmos(k))
          cfluxcggat(k)   = cfluxcgrow(ilmos(k),jlmos(k))
          cfluxcsgat(k)   = cfluxcsrow(ilmos(k),jlmos(k))
          canresgat(k)    = canresrow(ilmos(k),jlmos(k))
          xdiffusgat(k)   = xdiffus(ilmos(k))
          prbfrhucgat(k)  = prbfrhucgrd(ilmos(k))
          extnprobgat(k)  = extnprobgrd(ilmos(k))
          stdalngat(k)    = stdalngrd(ilmos(k))
          tcanoaccgat_m(k)= tcanoaccrow_m(ilmos(k),jlmos(k))
          uvaccgat_m(k)   = uvaccrow_m(ilmos(k),jlmos(k))
          vvaccgat_m(k)   = vvaccrow_m(ilmos(k),jlmos(k))
          vgbiomasgat(k)  = vgbiomasrow(ilmos(k),jlmos(k))
          gavgltmsgat(k)  = gavgltmsrow(ilmos(k),jlmos(k))
          gavglaigat(k)   = gavglairow(ilmos(k),jlmos(k))
          nppgat(k)       = npprow(ilmos(k),jlmos(k))
          nepgat(k)       = neprow(ilmos(k),jlmos(k))
          hetroresgat(k)  = hetroresrow(ilmos(k),jlmos(k))
          autoresgat(k)   = autoresrow(ilmos(k),jlmos(k))
          soilrespgat(k)  = soilresprow(ilmos(k),jlmos(k))
          rmgat(k)        = rmrow(ilmos(k),jlmos(k))
          rggat(k)        = rgrow(ilmos(k),jlmos(k))
          nbpgat(k)       = nbprow(ilmos(k),jlmos(k))
          litresgat(k)    = litresrow(ilmos(k),jlmos(k))
          socresgat(k)    = socresrow(ilmos(k),jlmos(k))
          gppgat(k)       = gpprow(ilmos(k),jlmos(k))
          dstcemlsgat(k)  = dstcemlsrow(ilmos(k),jlmos(k))
          litrfallgat(k)  = litrfallrow(ilmos(k),jlmos(k))
          humiftrsgat(k)  = humiftrsrow(ilmos(k),jlmos(k))
          rmlgat(k)       = rmlrow(ilmos(k),jlmos(k))
          rmsgat(k)       = rmsrow(ilmos(k),jlmos(k))
          rmrgat(k)       = rmrrow(ilmos(k),jlmos(k))
          burnfracgat(k)  = burnfracrow(ilmos(k),jlmos(k))
          probfiregat(k)  = probfirerow(ilmos(k),jlmos(k))
          lucemcomgat(k)  = lucemcomrow(ilmos(k),jlmos(k))
          lucltringat(k)  = lucltrinrow(ilmos(k),jlmos(k))
          lucsocingat(k)  = lucsocinrow(ilmos(k),jlmos(k))
          dstcemls3gat(k) = dstcemls3row(ilmos(k),jlmos(k))
          faregat(k)      = farerow(ilmos(k),jlmos(k))
          gavgscmsgat(k)  = gavgscmsrow(ilmos(k),jlmos(k))
c  
100   continue
c
      do 101 l=1,icc
       do 101 k=1,nml
          fcancmxgat(k,l)  = fcancmxrow(ilmos(k),jlmos(k),l)
          ailcggat(k,l)    = ailcgrow(ilmos(k),jlmos(k),l)
          ailcgsgat(k,l)   = ailcgsrow(ilmos(k),jlmos(k),l)
          fcancsgat(k,l)   = fcancsrow(ilmos(k),jlmos(k),l)
          fcancgat(k,l)    = fcancrow(ilmos(k),jlmos(k),l)
          co2i1cggat(k,l)  = co2i1cgrow(ilmos(k),jlmos(k),l)
          co2i1csgat(k,l)  = co2i1csrow(ilmos(k),jlmos(k),l)
          co2i2cggat(k,l)  = co2i2cgrow(ilmos(k),jlmos(k),l)
          co2i2csgat(k,l)  = co2i2csrow(ilmos(k),jlmos(k),l)
          slaigat(k,l)     = slairow(ilmos(k),jlmos(k),l)
          ancsveggat(k,l)  = ancsvegrow(ilmos(k),jlmos(k),l)
          ancgveggat(k,l)  = ancgvegrow(ilmos(k),jlmos(k),l)
          rmlcsveggat(k,l) = rmlcsvegrow(ilmos(k),jlmos(k),l)
          rmlcgveggat(k,l) = rmlcgvegrow(ilmos(k),jlmos(k),l)
          anveggat(k,l)    = anvegrow(ilmos(k),jlmos(k),l)
          rmlveggat(k,l)   = rmlvegrow(ilmos(k),jlmos(k),l)
          pfcancmxgat(k,l) = pfcancmxrow(ilmos(k),jlmos(k),l)
          nfcancmxgat(k,l) = nfcancmxrow(ilmos(k),jlmos(k),l)
          stemmassgat(k,l) = stemmassrow(ilmos(k),jlmos(k),l)
          rootmassgat(k,l) = rootmassrow(ilmos(k),jlmos(k),l)
          gleafmasgat(k,l) = gleafmasrow(ilmos(k),jlmos(k),l)
          bleafmasgat(k,l) = bleafmasrow(ilmos(k),jlmos(k),l)
          ailcbgat(k,l)    = ailcbrow(ilmos(k),jlmos(k),l)
          flhrlossgat(k,l) = flhrlossrow(ilmos(k),jlmos(k),l)
          pandaysgat(k,l)  = pandaysrow(ilmos(k),jlmos(k),l)
          lfstatusgat(k,l) = lfstatusrow(ilmos(k),jlmos(k),l)
          grwtheffgat(k,l) = grwtheffrow(ilmos(k),jlmos(k),l)
          lystmmasgat(k,l) = lystmmasrow(ilmos(k),jlmos(k),l)
          lyrotmasgat(k,l) = lyrotmasrow(ilmos(k),jlmos(k),l)
          tymaxlaigat(k,l) = tymaxlairow(ilmos(k),jlmos(k),l)
          stmhrlosgat(k,l) = stmhrlosrow(ilmos(k),jlmos(k),l)
          bmasveggat(k,l)  = bmasvegrow(ilmos(k),jlmos(k),l)
          rothrlosgat(k,l) = rothrlosrow(ilmos(k),jlmos(k),l)
          veghghtgat(k,l)  = veghghtrow(ilmos(k),jlmos(k),l)
          rootdpthgat(k,l) = rootdpthrow(ilmos(k),jlmos(k),l)
          tltrleafgat(k,l) = tltrleafrow(ilmos(k),jlmos(k),l)
          tltrstemgat(k,l) = tltrstemrow(ilmos(k),jlmos(k),l)
          tltrrootgat(k,l) = tltrrootrow(ilmos(k),jlmos(k),l)
          leaflitrgat(k,l) = leaflitrrow(ilmos(k),jlmos(k),l)
          roottempgat(k,l) = roottemprow(ilmos(k),jlmos(k),l)
          afrleafgat(k,l)  = afrleafrow(ilmos(k),jlmos(k),l)
          afrstemgat(k,l)  = afrstemrow(ilmos(k),jlmos(k),l)
          afrrootgat(k,l)  = afrrootrow(ilmos(k),jlmos(k),l)
          wtstatusgat(k,l) = wtstatusrow(ilmos(k),jlmos(k),l)
          ltstatusgat(k,l) = ltstatusrow(ilmos(k),jlmos(k),l)
          nppveggat(k,l)   = nppvegrow(ilmos(k),jlmos(k),l)
          rmlvegaccgat(k,l)= rmlvegaccrow(ilmos(k),jlmos(k),l)
          rmsveggat(k,l)   = rmsvegrow(ilmos(k),jlmos(k),l)
          rmrveggat(k,l)   = rmrvegrow(ilmos(k),jlmos(k),l)
          rgveggat(k,l)    = rgvegrow(ilmos(k),jlmos(k),l)
          gppveggat(k,l)   = gppvegrow(ilmos(k),jlmos(k),l)
          nepveggat(k,l)   = nepvegrow(ilmos(k),jlmos(k),l)
          nbpveggat(k,l)   = nbpvegrow(ilmos(k),jlmos(k),l)
          autoresveggat(k,l)= autoresvegrow(ilmos(k),jlmos(k),l)
          ailcmingat(k,l)  = ailcminrow(ilmos(k),jlmos(k),l)
          ailcmaxgat(k,l)  = ailcmaxrow(ilmos(k),jlmos(k),l)
          vgbiomas_veggat(k,l) = vgbiomas_vegrow(ilmos(k),jlmos(k),l)
          pftexistgat(k,l) = pftexistrow(ilmos(k),jlmos(k),l)

c         fire emission variables
          emit_co2gat(k,l)  = emit_co2row(ilmos(k),jlmos(k),l)  
          emit_cogat(k,l)   = emit_corow(ilmos(k),jlmos(k),l)   
          emit_ch4gat(k,l)  = emit_ch4row(ilmos(k),jlmos(k),l)  
          emit_nmhcgat(k,l) = emit_nmhcrow(ilmos(k),jlmos(k),l) 
          emit_h2gat(k,l)   = emit_h2row(ilmos(k),jlmos(k),l)   
          emit_noxgat(k,l)  = emit_noxrow(ilmos(k),jlmos(k),l)  
          emit_n2ogat(k,l)  = emit_n2orow(ilmos(k),jlmos(k),l)  
          emit_pm25gat(k,l) = emit_pm25row(ilmos(k),jlmos(k),l) 
          emit_tpmgat(k,l)  = emit_tpmrow(ilmos(k),jlmos(k),l)  
          emit_tcgat(k,l)   = emit_tcrow(ilmos(k),jlmos(k),l)   
          emit_ocgat(k,l)   = emit_ocrow(ilmos(k),jlmos(k),l)   
          emit_bcgat(k,l)   = emit_bcrow(ilmos(k),jlmos(k),l) 
          burnveggat(k,l)   = burnvegrow(ilmos(k),jlmos(k),l) 
c
101   continue
c
      do 102 l=1,iccp1
       do 102 k=1,nml
          litrmassgat(k,l)=litrmassrow(ilmos(k),jlmos(k),l)
          soilcmasgat(k,l)=soilcmasrow(ilmos(k),jlmos(k),l)
          hetroresveggat(k,l)= hetroresvegrow(ilmos(k),jlmos(k),l)
          litresveggat(k,l)= litresvegrow(ilmos(k),jlmos(k),l)
          soilcresveggat(k,l)= soilcresvegrow(ilmos(k),jlmos(k),l)

102   continue
c
      do 105 l=1,12     !12 months
       do 105 k=1,nml
          mlightnggat(k,l)=mlightnggrd(ilmos(k),l)
105   continue
c
      do 106 l=1,2     !2 pfts (ndl dcd & crops)
       do 106 k=1,nml
          colddaysgat(k,l)=colddaysrow(ilmos(k),jlmos(k),l)
106   continue
c
      do 201 l=1,ican
       do 201 k=1,nml
          ailcgat(k,l)     = ailcrow(ilmos(k),jlmos(k),l)
          zolncgat(k,l)    = zolncrow(ilmos(k),jlmos(k),l)
          cmasvegcgat(k,l) = cmasvegcrow(ilmos(k),jlmos(k),l)
          paicgat(k,l)     = paicrow(ilmos(k),jlmos(k),l)
          slaicgat(k,l)    = slaicrow(ilmos(k),jlmos(k),l)
c
          alvsctmgat(k,l)  = alvsctmrow(ilmos(k),jlmos(k),l)
          alirctmgat(k,l)  = alirctmrow(ilmos(k),jlmos(k),l)
201   continue
c
      do 250 l=1,ignd
       do 250 k=1,nml
          sandgat(k,l) = sandrow(ilmos(k),jlmos(k),l)
          claygat(k,l) = clayrow(ilmos(k),jlmos(k),l)
          orgmgat(k,l) = orgmrow(ilmos(k),jlmos(k),l)
          tbaraccgat_m(k,l) = tbaraccrow_m(ilmos(k),jlmos(k),l)
250   continue
c
      do 280 l=1,icc
       do 280 m=1,ignd
        do 280 k=1,nml
          rmatctemgat(k,l,m) = rmatctemrow(ilmos(k),jlmos(k),l,m)
280   continue
c
      do 290 l=1,ican
       do 290 m=1,ignd
        do 290 k=1,nml
          rmatcgat(k,l,m) = rmatcrow(ilmos(k),jlmos(k),l,m)
290   continue
c

      return
      end
