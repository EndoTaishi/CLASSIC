	 subroutine mosspht(ilg,ignd,isand,iday,qswnv,thliq,tbar,thpor,
	1		co2conc,tsurfk,zsnow,delzw,pres,qg,coszs,Cmossmas,dmoss,
c		output below
	2		anmoss,rmlmoss,cevapms,ievapms,ipeatland	
c	testing
	3		,iyear, ihour,imin ) 
C	
c	-----------definitions---------------------------------------------
	 implicit none

	 integer i, j, ilg,ignd, isand(ilg,ignd), iday, ihour, iyear,imin
	 integer ievapms(ilg),        ipeatland(ilg)
	 real:: qswnv(ilg)	!visable short wave radiation = qswnv in TSOLVE
	                    !and qswnvg in TSOLVC (W/m2)
c	environment related
	 real	thliq(ilg,ignd), 	tbar(ilg,ignd),	thpor(ilg,ignd),	
	1		bi(ignd),			co2conc(ilg), 		
	2		zsnow(ilg),		delzw(ilg,ignd),	pres(ilg),
	3		fcs(ilg),    	fgs(ilg),    	fc(ilg),     	fg(ilg),
	4		qg(ilg),		coszs(ilg)
c	moss related
	 real:: Cmossmas(ilg)	!unig kg moss C updated in ctem 
	 real:: mmoss(ilg)		!unit kg moss mass in kg, update later 	
	 real:: dmoss(ilg)		!unit m, depth of living moss. assume = 2 cm
						!can be related to mmoss as a variable	
c	---------------input above this line-------------------------------

	 real:: anmoss(ilg)		!net photosynthesis (umol/m2/s)	
	 real:: rmlmoss(ilg)	!moss autotrophic respiration (umol/m2/s)

c	---------------output above this line, parameters below------------

	 real:: rmlmoss25 = 1.1	!base dark respiration rate of moss umol/m2/s 
						!(best estimate is from fig. 6b in 
						!Williams and Flanagan (1998)
						!seasonal variation not considered		
	 real:: tau25m=2321.0	!tau coefficient (rate at 25 Celcius) is 
						!computed on the basis of the specificity 
						!factor (102.33) times Kco2/Kh2o (28.38) 
						!to convert for value in solution to that 
						!based in air. The old value was 2321.1.
						!New value (2904.12) Quercus robor (Balaguer
						! et al., 1996) !Similar number from Dreyer 
						!et al. 2001, Tree Physiol, tau= 2710 
						!2600 from Brooks and Farquhar (1985) used 
						!AS A CONSTANT; 2321 Harley et al. (1992) 
	 real:: tref = 298.16	!unit K
	 real:: ektau =-29000.0	!J mol-1 (Jordan and Ogren, 1984)
	 real:: gasc = 8.314	!gas constant (J mol-1 K-1
	 real:: kc25= 30.2 		!kinetic coef for CO2 at 25 C (Pa)
						!30.2 Arbutus unedo  (Harley et al., 1986)
						!40.4 (von Caemmerer et al., 1994)
						!27.46 deciduous forest understory (Wilson
						! et al. 2000); 30.0 (Collatz et al., 1991)
						!30.0 Populus tremula,Coylus avelana,Tilia
						!cordata (Kull and Niinemets, 1998)
	 real:: ko25=24800.0	!kinetic coef for O2 at 25C (rate)  (Pa)
						!24800 Arbutus unedo (Harley et al.,1986)
						!41980 deciduous forest understory (Wilson et al. 2000)
						!30000 (Collatz et al., 1991)
						!30000 (Kull and Niinemets, 1998) see kc25
      real:: ec=63500.0		!Activation energy for Kc of CO2 (J mol-1) at 25C
						!59430 Arbutus unedo (Harley et al., 1986) 
						!80500 deciduous forest understory (Wilson et al., 2000)
						!63500 mosses and shrubs Tenhunen et al. (1994)
	 real:: ej=37000.0		!activation energy for electron transport, 
						!&37000 at 25C J mol-1(Farquhar et al.(1980)
						!55000 J mol-1 at 311K  (Harley & Baldocchi, 
						!1995,PCE) lichis tree
      real:: eo=35000.0		!activation energy for Ko (J/mol) at 25C
						!36000 Arbutus unedo (Harley et al.,1986)
						!14500 for deciduous forest understory (Wilson et al., 2000)
						!35000 mosses and shrubs (Tenhunen et al.,1994) 
	 real:: evc=53000.0		!activation energy for carboxylation at 25C, J mol-1
						!53000 (Kirschbaum & Farquhar, 1984)
						!55000 (at 311K) (Harley & Baldocchi, 1995,PCE)  

	 real:: sj=719.0		!constant affecting J at low temperature  
						!at 25C(j/mol); 714.0 in Lloyd et al. (1995)  
						!Macadamia integrifolia and Litchi chinesis
	 real:: hj=220300.0		!constant affecting J at high temperature 
						!at 25C(j/mol); 220300 Macadamia integrifolia 
						!and Litchi chinesis Lloyd et al. (1995)
						!206083 Pinus radiata (Wilson et al. 2000)
	 real:: alpha=0.21		!efficiency of conversion of incident photons
						!into electrons(mol electron/mol photon), 
						!not sensitive, alpha = 0.21 Flanagan's 
						!value for sphagnum and pleurozium
						!!alpha = 0.3 (Seel et al.,1992) for
						!T. ruraliformis and d. palustris
	 real:: std_press = 101325.0	!standard atmospheric pressure 

c
c	constants / parameters above, temporary variables below--------------
c
	 integer:: pheno(ilg) 	!phenology flag of mosses, 1 = photosynthesis, 
						!0 = does not photosynthesis 
	 real:: parm(ilg) 		!par at the ground (moss layer) surface umol/m2/s
	 real:: tsurf(ilg)		!grid average ground surface temperature in C
	 real:: tsurfk(Ilg) 	!grid average ground surface temprature in K
	 real:: wmoss(ilg) 		!water content extraporated from the surface
						!humidity qg and thliq of the first soil layer
						!unit g water/ g dw 
	 real:: wmosmin(ilg)	!residual water content kg water /kg moss 
	 real:: wmosmax(ilg)	!maximum water content kg water /kg moss 
	 real:: fwmoss(ilg)		!relative water content of mosses in g fw /g dw
	 real:: dsmoss(ilg)		!degree of moss saturation = relative water 
						!content/maximum relative water content 
	 real:: cevapms(ilg)	!evaporation coefficent for moss surface	 
	 real:: g_moss(ilg)		!moss conductance umol/m2/s (based on 
						!Williams and Flanagan,1998 for Sphagnum)
	 real:: mwce (ilg)		!moisture function of dark respiration of moss
	 real:: tmoss(ilg)		!moss temperature extraporated from the tbar 1
						!and grid averaged ground surface temperature
						!tsurf
	 real:: tmossk(ilg)		!moss temperature in K
	 real:: q10rmlmos(ilg)	!temperature function of the moss dark respiration 
	 real:: gamma(ilg)		!compensation point for gross photosynthesis (Pa)	
	 real:: o2(ilg)		!partial presure of oxygen (Pa)
	 real:: co2a(ilg)		!partical pressure of co2 (pa) same as in PHTSYN
	 real:: tau(ilg)		!arrhenius funciton of temperature
	 real:: kc(ilg)		!kinetic coeffficient of CO2 for photosynthesis
	 real:: ko(ilg)		!				.....o2... 
	 real:: bc(ilg)		!coefficient used for wc
	 real:: vcmax25(ilg)	!seasonal varied maximum carboylation at 25  
						!sphagnum (fig. 6, Williams and Flanagan, 1998)
	 real:: vcmax(ilg)		!max carboxylation rate (umol/m2/s)
	 real:: jmax25(ilg)		!maximum electorn transport rate at 25 degrees (umol/m2/s)
	 real:: jmax(ilg)		!maximum electorn transport rate (umol/m2/s)
	 real:: wj(ilg)		!net co2 assimilation rate limited by 
						!electron transport (umol/m2/s)=jE in PHTSYN
	 real:: wc(ilg)		!net co2 assimilation rate limited by 
	 real:: ws(ilg)		!net co2 assimilation rate limited by 
						!sucrose availability (umol/m2/s)=JE IN PHTSYN
	 real:: photon(ilg)		!electron transport rate (umol/m2/s)
c	--------------internal variables above this line------------------

	 real:: term1(ilg), term2(ilg), term3(ilg)	!temporal terms for 
						!photosynthesis calculations
	 real:: psna(ilg), psnb(ilg), psne(ilg)	!coefficients for quadratric 
						!solution of net photosynthesis
	 real:: mI(ilg), mII(ilg)!coefficients of the solutions for net psn
	 real:: beta1, beta2
	 real:: temp_b, temp_c, temp_r, temp_q1, temp_q2, temp_jp
c	 -----------temporal terms above this line------------------------

	 real     DELT,TFREZ, HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,SPHW,
	1	     SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP
	
	 real     zolnms,thpms,thrms,thmms,bms,psisms,grksms,hcpms,
	1		sphms,rhoms,slams

c	-------------common block parameters above this line--------------
 
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
	 common /peatland/ zolnms,thpms,thrms,thmms,bms,psisms,grksms,
	1				hcpms, sphms,rhoms,slams	
C     PHOTOSYNTHESIS COUPLING OR CURVATURE COEFFICIENTS
	 DATA BETA1/0.950/
      DATA BETA2/0.990/
          
c	find the light level (parm) at the ground surface for moss photosynthesis,  
c    and a scaling factor degree of saturation of the moss layer phenology
c	parm is in umol/m2/s and converted from qswnv in W/m2

	 do 100 	i = 1, ilg
		parm(i)= qswnv(i)*4.6
		o2(i)  = 20.9/100.0 * pres(i)
		co2a(i) = co2conc(i)/1000000.0 * pres(i)		 
		wj(i)	= 0.0
		ws(i)	= 0.0
		wc(i)	= 0.0
		anmoss(i) = 0.0
		rmlmoss(i)= 0.0
		tsurf(i)= tsurfk(i)-tfrez    		
		tmossk(i) = tsurfk(i)
		tmoss(i)= tsurf(i)
100	 continue
c
c	phenology  water factor on mosses ,grow when temperature > -4.0 
c	and snowpack < 0.15 m
	 do 200 i = 1, ilg
		if (zsnow(i) .gt. 0.15 .or. tsurf (i) .lt. -4.0)  then		
			pheno(i) = 0.0
		else
			pheno(i) = 1.0
		endif	
200	 continue
c
c    
c	** water content used for the living moss (depth dmoss)
c    dmoss is an input and site specific. Preferably make dmoss a function  
c    of Cmoss and ipeatland (different species in fens and bogs)
c 	observed range of wmoss: 5 to 40 in Robrek (2007, 2009), 5 to 25 
c    (Flanagen and Williams 1998)

	 do 250 	i = 1, ilg
		mmoss(i) = Cmossmas(i)/0.46
		wmoss(i)= thliq(i,1)*rhow/(mmoss(i)/dmoss(i))
		wmosmax(i) = min(45.0, thpms*dmoss(1)*rhow/mmoss(i)) 
		wmosmin(i) = max(5.0, thmms*dmoss(1)*rhow/mmoss(i)) 
		wmoss(i)= min(wmosmax(i),max(wmosmin(i),wmoss(i)))
		fwmoss(i)=wmoss(i)+1	!g fresh weight /g dry weight
250	 continue
c
c	** moss conductance g_moss umol/m2/s 
c    (Williams and Flanagan, 1998 for Sphagnum). follow MWM, fwmoss is
c     the mosswat_fd in MWM. Empirical equation is only valid up to 
c    fwmoss=13, above 13 apply a linear extension to the equation.

	 do 300 i = 1, ilg
		if (fwmoss(i) .le. 13.0) 	then
	  	  g_moss(i)=-0.195+0.134*fwmoss(i) - 0.0256*(fwmoss(i))
     1		**2 + 0.00228*(fwmoss(i))**3 - 0.0000984*
     2		(fwmoss(i))**4 + 0.00000168*(fwmoss(i))**5
	  	else	
		  g_moss (i) = -0.000447 * fwmoss(i) + 0.0489	
		endif				
		  g_moss (i) = g_moss(i) * 1000000.0
		  g_moss(i) = max(0.0, g_moss(i)) 

c	** moss surface evaporation coefficient 
c    controled by the degree of saturation in moss, pass to TSOLVC and TSOLVE
C	apply a similar equation of soil surface cevap in TPREP 
c    CEVAP = 0.25*[1 – cos(THLIQ*pi/THFC)]^2

		 dsmoss(i) = (wmoss(i)-wmosmin(i))/(wmosmax(i)-wmosmin(i))
		 if (dsmoss(i) .lt.   0.001)       then
			ievapms(i) = 0
			cevapms(i) = 0.
		 elseif (dsmoss(i) .ge. 1.0)		then
			ievapms(i) = 1
			cevapms(i) = 1.0	
		 else 
			ievapms(i) = 1
			cevapms(i) = 0.25*(1.0-cos(3.14159*dsmoss(i)))**2
		 endif

c
c	** moss water content effect on dark respiration 
c	in MWM and PDM an optimal wmoss is at 5.8 gw/gdw(fig. 2e, Frolking et al.,1996)
c    Recent studies show weak but significant increases of sphagnum dark respiration
c    with moss water content above 5.8 gw/gdw (Adkinson and Humphreys, 2011 and ref.)
c    this change has improved the ER simulation greatly
		if (wmoss(i) .lt. 0.4)	then
			mwce (i) = 0.0
	 	elseif(wmoss(i) .lt. 5.8 .and. wmoss(i) .gt. 0.4) then
			mwce(i) = 0.35*wmoss(i)**(2.0/3.0)-0.14
  		else
c			mwce(i) =-0.04*wmoss(i)+1.232 ! in MWM and PDM
               mwce(i)= 0.01*wmoss(i)+0.942 
		endif
300	 continue
c
c	** moss dark respiration
c     observed range of rmlmoss 0.60 to 1.60 umol/m2/s (e.g. Adkinson 2006)
	 do 350 i = 1, ilg
	     q10rmlmos(i)=(3.22-(0.046*tmoss(i)))**((tmoss(i)-25.0)/10.0)
		rmlmoss(i) = rmlmoss25*mwce(i)*q10rmlmos(i)

c	** moss photosynthesis
c    calculate bc (coefficient used for Wc, limited by Rubisco)
c
		tau(i) = tau25m*exp((tmossk(i)-tref)*ektau/
     1		     (tref*gasc*tmossk(i)))
		gamma(i) = .5 * o2(i)/ tau(i)		
		kc(i) = kc25*exp((tmossk(i)-tref)*ec/(tref*gasc*tmossk(i)))
		ko(i) = ko25*exp((tmossk(i)-tref)*eo/(tref*gasc*tmossk(i)))
		bc(i)  = kc(i) * (1.0 + (o2(i)/ ko(i)))	
c		
c    seasonal change of Vcmax, sphagnum (fig. 6, Williams and Flanagan, 1998)
c
		if(iday .lt. 121)        then
			vcmax25(i) = 6.0
		else if(iday .gt. 245)   then
			vcmax25(i) = 7.0
		else 
			vcmax25(i) = 13.5 
		endif	
		vcmax(i) = vcmax25(i)*exp((tmossk(i)-tref)*evc/
	1			(tref*gasc*tmossk(i)))
c
C     calculate ws (phototysnthesis rate limited by transport capacity)
C	 = js in PHTSYN3
c
          if (coszs(i).gt.0.0)	THEN
			ws(i) = 0.5*vcmax(i)
		endif
c
c	calculate the maximum electorn transport rate Jmax (umol/m2/s)
c    1.67 = vcmax25m/jmax25m ratio
c
		jmax25(i) = 1.67 * vcmax25(1)	 
		term1(i)=exp(((tmossk(i)/tref)-1)*ej/(gasc * tmossk(i)))
     	term2(i)=1+exp(((tref*sj)-hj)/(tref*gasc))
		term3(i)=1+exp(((sj*tmossk(i))-hj)/(gasc*tmossk(i)))
		jmax(i)=jmax25(i)*term1(i)*term2(i)*term3(i)
		if (jmax(i) .gt. 0.0) 			then
			photon(i)=alpha*parm(i)/sqrt(1.0+(alpha**2*parm(i)**2		 
     &	     		/(jmax(i)**2)))	
		else
			photon(i) = 0.0	!electron trasport rate in mosses 					
		end if
c
c	calculate Wj, Wc (Farquhar and Caemmerer 1982)
c	wj = light limited, = je in PHTSYN
c
		wj(i) = photon(i)*(co2a(i)-gamma(i))/
	1			(4.*co2a(i)+(8*gamma(i)))	
c
c	carboxylase(rubisco) limitation = jc in PHTSYN
		wc(i) = vcmax(i)*(co2a(i)-gamma(i))/(co2a(i)+bc(i))		

c	Choose the minimum of Wj and Wj both having the form: 
c	W = (a Ci - ad) / (e Ci + b)
c	Then set a, b, d and e for the quadratic solution for net photosynthesis. 
c
		if(wj(i) < wc(i))then
			psnb(i) = 8 * gamma(i)
			psna(i) = photon(i)
			psne(i) = 4.0
		elseif(wc(i) < wj(i))then 		
			psnb(i) = bc(i)  
			psna(i) = vcmax(i)
			psne(i) = 1.0
		end if
350    continue
c
c	Calculate net and gross photosynthesis by solve the quadratic equation
c	first root of solution is net photosynthesis An= min(Wj,Wc) - Rd 
c	gross photosynthesis GPP = min(Wc,Wj) = An + Rd 
	
    	 do 400	i = 1, ilg
		if (psna(i) .gt.0.0) 			               then 		
		    mI(i)=rmlmoss(i)-(psnb(i)*g_moss(i)/pres(i)/psne(i)) 
     &			- (co2a(i)*g_moss(i)/pres(i))-(psna(i)/psne(i))	
		    mII(i)=(psna(i)*co2a(i)*g_moss(i)/pres(i)/psne(i))-
  	1		(rmlmoss(i)*co2a(i)*g_moss(i)/pres(i))-(rmlmoss(i)*
   	2 		psnb(i)*g_moss(i)/pres(i)/psne(i))-(psna(i)*gamma(i)*
	3		g_moss(i)/pres(i)/psne(i))	
		else 	
     	    mI(i)= 0.0
		    mII(i)= 0.0
		endif		
		anmoss(i) = (-mI(i)-(mI(i)*mI(i)-4*mII(i))**0.5)/2	
		anmoss(i) = min(anmoss(i), ws(i))		
400		continue

c    write to output file midday             CT16D_G 
	 if (iyear .eq. 2004 .and. ihour==12)        then
	  	write(98,6998)	iday,tmoss,cevapms,fwmoss,thliq(1,1), 
     1    dsmoss,g_moss, wmoss,rmlmoss,mwce,q10rmlmos,wmosmax,wmosmin
     
6998   	format(I6,15f12.4) 		
	 endif

	 return
	 end


