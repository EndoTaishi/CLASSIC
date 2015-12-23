	  subroutine  decp(il1,il2,iyear,iday,ihour,imin,ipeatland,isand,
	1	litrmassms,hpd, wtable,tbar, thliq, thice,thpor,bi,zbotw,
	2    delzw,psisat,tfrez, jdst,jdend,iyd,
c    -------------- inputs above this line, outputs below -------------
     3    litresms, socresp, resoxic, resanoxic)  

c     grid average peat soil heterotrophic respiration subroutine 
c                        YW March 20, 2015 

      use ctem_params,      only :icc, ilg,ignd,zero,tanhq10

	 implicit none
c	 inputs-----------------------------------------------------------
 	 integer  iyear, iday, ihour, imin ,i,j, il1, il2, isand(ilg,ignd)
 	 integer  jdst, jdend, iyd
	 integer	ipeatland(ilg)	!0 = not peatland, 1 = bog, 2 = fen  
      integer:: lewtable(ilg) !layer index of the water table layer 
	 real     hpd(ilg) , wtable(ilg), tbar(ilg,ignd),	!(K) 		 
	1         thliq(ilg,ignd),	thice(ilg,ignd),    thpor(ilg,ignd), 
	2	     bi(ilg,ignd),       zbotw(ilg,ignd),	delzw(ilg,ignd),	
	3         frac(ilg),          tfrez,              litrmassms(ilg)	
c	 --------------outputs C fluxes in umolCO2/m2/s--------------------
 	 real:: litresms(ilg)   !moss litter respiration 
	 real:: socresp(ilg)    !soil C respiration 
	 real:: resoxic(ilg)    !respiration rate of the oxic compartment
	 real:: resanoxic(ilg)  !respiration rate of the oxic compartment
c	 internal variables------------------------------------------------ 
	 real:: Cso (ilg) 	!carbon mass in the oxic compartment (kgC/m2)
	 real:: Csa (ilg)	!carbon mass in the anxic compartment (kgC/m2) 
	 real:: fto(ilg) 	!temperature factor of the oxic soil respiration
	 real:: fta(ilg) 	!temperature factor of the anoxic soil respiration
	 real:: tsoilo(ilg)	!average temperature in the oxic compartment(C)
	 real:: tsoila(ilg) !average temperature of the anoxic compartment (C)
	 real:: ewtable(ilg)!effective water table depth (m)
	 real:: ratescpo(ilg) !oxic respiration rate constant(umolCO2/kgC/s)
	 real:: ratescpa(ilg) !anoxic respiration rate constant(umolCO2/kgC/s)
	 real:: psi(ilg,ignd) !matrix potential of soil layers (Pa)	
	 real:: psisat(ilg,ignd) !saturated matrix potential in soil (m)
	 real:: q10funcms(ilg)	!q10 fuction for moss litter respiration
	 real:: ltrmosclms(ilg)	!moisture scale for litter respiration
	 real:: litrtempms(ilg)	!temperature of the moss litter
	 real:: litpsims(ilg) 	!matrix potential of litter 	 
	 real:: litrq10ms(ilg)	!q10 coefficient as a functon of T (as in CTEM)
	 real:: soilq10o(ilg)	!q10 coefficient of oxic soil
	 real:: soilq10a(ilg)	!q10 coefficient of anoxic soil
							
c	parameters (from PCARs) can be moved to ctem_params----------------
	 real::	dctmin= 269.0	!minimum temperature of soil respiration,K
	 real:: 	dcbaset=283.0	!base temperature for Q10, K 
	 real:: 	bsrateltms = 0.20	!yr-1, 0.05 (T.Moore unpubished data)
							!for shrubs 0.2 from Moore and around
							!0,5 for ctem. TM's values is natural
							!condition, CTEM is opt at 15degrees,
							!So TM's values can be scaled up
							!more rates in Aerts 1997, Moore 2007
c	------------------------------------------------------------------
c
c	initialization
	 do 10 i = il1, il2
		litresms(i) = 0.0
		socresp(i) = 0.0 
		resoxic(i) = 0.0	
		resanoxic(i) = 0.0	
10	 continue
c
c    ** calculate soil respiration from peat 
c
c    find the effective water table depth and the layer index to devide
c    the peat soil into two compartment 
c
	 do 20 	i= il1, il2
          ewtable(i)= wtable(i)            
          if (ewtable(i) .le. 0.0)                then
              lewtable(i) = 0          
          elseif (ewtable(i) .le. 0.1)            then
              lewtable(i) = 1
          else        
              do j = 1, ignd
                  if (ewtable(i) .gt. zbotw(i,j))        then
                     lewtable(i)=j+1
                  endif  
              enddo
          endif
c              
c    find the temperature in litter, oxic soil and anoxic soil in kalvin
c    lewtable is the layer index of the water table layer, lewtable = 0 
c    indicates WTD is above the ground surface. 
c    Set the oxic layer temperature to dctmin (minimum soil respiration
c    temperture) when the entire soil is in the anoxic zone. 
c
          tsoilo(i) = 0.0
          tsoila(i) = 0.0           
          if (lewtable(i) .eq. 0)   then  !WT is at or above the surface
               do j = 1, ignd
                   tsoila(i) = tsoila(i)+tbar(i,j)*delzw(i,j)
               enddo
               tsoila(i)=tsoila(i)/zbotw(i,ignd)
               tsoilo(i) = dctmin	
          elseif (lewtable(i) .eq. 1) then	!WT is at the first layer
               tsoilo(i)=tbar(i,1)
               do j= lewtable(i)+1, ignd
                   tsoila(i)=tsoila(i)+ tbar(i,j)*delzw(i,j)
               enddo
               tsoila(i)=(tsoila(i)+tbar(i,1)*(zbotw(i,lewtable(i))
     1             -ewtable(i)))/(zbotw(i,ignd)-ewtable(i))
          else	                             !WT is below layer 1
              do j = 1,lewtable(i)-1
                  tsoilo(i) = tsoilo(i)+ tbar(i,j)*delzw(i,j)
              enddo
              tsoilo(i)= (tsoilo(i)+ tbar(i,lewtable(i))*
     1               (ewtable(i)-zbotw(i,lewtable(i)-1)))/ewtable(i)
              do j = ignd,lewtable(i)+1,-1
                  tsoila(i) = tsoila(i)+ tbar(i,j)*delzw(i,j)
              enddo
              tsoila(i)= (tsoila(i)+tbar(i,lewtable(i))*
     1               (zbotw(i,lewtable(i))-ewtable(i)) )/
     2               (zbotw(i,ignd)-ewtable(i))
          endif
20	 continue
c
c	calculate the temperature multiplier (ftsocres) for oxic and anoxic
c    soil compartments
c
	 do 30 	i = il1, il2	
        	soilq10o(i) = tanhq10(1) + tanhq10(2)*
     &            ( tanh( tanhq10(3)*(tanhq10(4)-(tsoilo(i)-tfrez))))
		fto(i)= soilq10o(i)**(0.1*(tsoilo(i)-tfrez-15.0))
        	soilq10a(i) = tanhq10(1) + tanhq10(2)*
     &            ( tanh( tanhq10(3)*(tanhq10(4)-(tsoila(i)-tfrez))))
		fta(i)= soilq10a(i)**(0.1*(tsoila(i)-tfrez-15.0))
30	 continue
c
c	find the heterotrophic respiration rate constant in tje oxic and
c    anoxic (unit in yr-1), based on Fig.2b in Frolking 2001 
c		
	 do 40 i = il1, il2
	     if (ipeatland(i) == 1) 			then		!bogs
	   	    if (ewtable(i) .lt. 0.0) 		                    then
		        ratescpo(i)=0.0
		        ratescpa(i)=-0.183*exp(-18.0*hpd(i))+0.03
	1			          *hpd(i)+0.0134
		    elseif (ewtable(i).lt.0.30 .and.ewtable(i).ge.0.0) then
		        ratescpo(i)=0.009*(1-exp(-20*ewtable(i)))
     1                        +0.015*ewtable(i)
		        ratescpa(i)=0.009*exp(-20*ewtable(i))-0.183*exp(-18*
	1			          hpd(i))-0.015*ewtable(i)+0.0044
	   	    elseif (ewtable(i) .ge. 0.30) 	                    then
         		   ratescpo(i)=0.0134-0.183*exp(-18*ewtable(i))+0.003*
	1			          ewtable(i)
		        ratescpa(i)=-0.183*exp(-18*hpd(i))+0.003*(hpd(i)
	1			          -wtable(i))+0.183*exp(-18*ewtable(i)) 
c          	   ratescpa(i)=-0.183*exp(-18*hpd(i))+0.003*(hpd(i)
c	1			   -wtable(i))+0.183*exp(-18*ewtable(i))-0.004504   !for continuity
	         endif 	 			     
	     elseif (ipeatland(i) == 2)			then		!fens	     
	  	    if (ewtable(i) .lt. 0.0) 		                    then		
		        ratescpo(i) = 0.0
		        ratescpa(i) = 0.01512 -1.12*exp(-25*hpd(i))   
	   	    elseif(ewtable(i).lt. 0.30 .and.ewtable(i).ge. 0.0)then
		        ratescpo(i) = -0.01*exp(-40*ewtable(i))+0.015*
	1				     ewtable(i)+0.01
		        ratescpa(i) = abs(-0.01*exp(-40*ewtable(i))-1.12*exp(
	1			          -25*hpd(i))+0.015*ewtable(i)+0.005119)
	   	    elseif(ewtable(i) .ge. 0.30) 	                    then
     		   ratescpo(i) = 0.01512-1.12*exp(-25*ewtable(i))
     		   ratescpa(i) = -1.12*(exp(-25*hpd(i))-exp
	1			          (-25*ewtable(i)))
	   	    endif  				
	     endif					
c	converts respiration rates from kg c/kg c.year to u-mol co2/kgC/s
		ratescpo(i) = 2.64 * ratescpo(i)
		ratescpa(i) = 2.64 * ratescpa(i) 	
	
c	find the carbon storage in oxic and anoxic compartments (Cso. Csa)
c	The water table depth delineates the oxic and anoxic compartments.
c	functions (R**2 = 0.9999) determines the carbon content of each 
c	compartment from a peat bulk density profile based on unpulished 
c	data from P.J.H. Richard	(described in fig. 1, Frokling et al.(2001)
c	conversion of peat into carbon with 48.7% (Mer Bleue unpublished data, 
c	Moore)

	     Cso(i) = (4056.6*ewtable(i)**2+72067.0*ewtable(i))
	1              *0.487/1000.0
	     Csa(i) = ((4056.6*hpd(i)**2+72067.0*hpd(i))*0.487/1000.0)
	1			-Cso(i)	

40	 continue

c	find the soil respiration rate in Cso and Csa umol/m2/s.
c	Moisture multiplier (0.025) indicates rate reduction in decomposition due 
c	to anoxia (Frolking et al. 2001), only applied to anoxic layer
	 do 50 i = il1, il2
		resoxic(i)   = ratescpo(i)*Cso(i)*fto(i)
		resanoxic(i) = ratescpa(i)*Csa(i)*fta(i) *0.025
		socresp(i)   = resoxic(i) + resanoxic(i)	
50	 continue

c	**calcualte litter respiration of moss

c    first find the matrix potential of the soil layers
	 do 60 j = 1, ignd
		do 60 i = il1, il2
             if(isand(i,j).eq.-3.or.isand(i,j).eq.-4) then  !ice or rock
	      	psi(i,j) = 10000.0 ! a large number so that ltrmoscl = 0.2
             else
               if (thliq(i,j)+ thice(i,j)+0.01 < thpor(i,j)
     1              .and.  tbar(i,j) <273.16)                   then
                  psi(i,j) = 0.001
               elseif (thice(i,j) > thpor(i,j))    then
                  psi(i,j) = 0.001   !set to saturation 
               else 
                  psi(i,j)=psisat(i,j)*(thliq(i,j)/(thpor(i,j)
     1                   -thice(i,j)))**(-bi(i,j))
                endif                                      
		   endif 
60	 continue	 
				        
c    litter in peatlands can be saturated so we limit the rate by high
c	moisuture level similar to soil in CTEM, but less effectively (the
c    min moisture factor is at 0.5 for moss litter but at 0.2 for soil). 
      do 70 i = il1, il2
		litpsims(i) = psi(i,1)
c    limit of ltrmoscalms at saturation YW April 10, 2015 
		if (litpsims(i) .gt. 10000.0) then
	    		ltrmosclms(i) = 0.2
        	elseif (litpsims(i).le. 10000.0 .and.litpsims(i).gt. 6.0) then
	          ltrmosclms(i)=1.0-0.8*((log10(litpsims(i))-log10(6.0))
	1				/(log10(10000.0)-log10(6.0)) )**1.
         	elseif (litpsims(i).le. 6.0 .and.litpsims(i).gt. 4.0) then
              	ltrmosclms(i)=1.0
        	elseif (litpsims(i).le. 4.0 .and. litpsims(i).gt.psisat(i,1)) 
	1												then 
              	ltrmosclms(i)=1.0-0.99*((log10(4.0)-log10(litpsims(i)))/   
     1         		(log10(4.0)-log10(psisat(i,1))))
         	elseif (litpsims(i) .le. psisat(i,1)) 				then
             	ltrmosclms(i)=0.01
        	endif
c    -----------------------------------------------------------------
c
         	ltrmosclms(i)=max(0.0,min(ltrmosclms(i),1.0))

c	find the temperature factor for moss litter respiration   
          litrtempms(i)=tbar(i,1)-tfrez	
        	litrq10ms(i) = tanhq10(1) + tanhq10(2)*
     &            ( tanh( tanhq10(3)*(tanhq10(4)-litrtempms(i))  ) )
		q10funcms(i)= litrq10ms(i)**(0.1*(litrtempms(i)-15.0))

c	calculate the litter respiration rate in mosses and converts it 
c     from kg c/kg c.year to u-mol co2/kg c.s using 2.64
       	litresms(i)=ltrmosclms(i)*litrmassms(i)*bsrateltms*
	1			2.64*q10funcms(i)
70		continue
c
c     write peat soil respiration details .CT15D_G
      
      if ((iyd.ge.jdst).and.(iyd.le.jdend)) then   

      write(97, 6997) litresms, litpsims(1), psisat(1,1),ltrmosclms, 
     1     litrmassms, tbar(1,1), q10funcms ,litrtempms,
     2     ratescpo, ratescpa, Cso,Csa, fto,fta,
     3     resoxic, resanoxic, frac, 
     4     tsoila-tfrez, tsoilo-tfrez,ewtable,real(lewtable),
     5     tbar(1,1)-tfrez, tbar(1,2)-tfrez, tbar(1,3)-tfrez,
     6     thliq(1,1) 
6997  format(50f10.3) 
      endif
c
 	 return
	 end 

