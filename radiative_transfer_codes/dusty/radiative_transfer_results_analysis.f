c     =========================================================================     
c     These are the subroutines analyzing the results of the radiative        
c     transfer calculation.                                   [MN, Mar'99]      
c     =========================================================================     
C     Table of Contents                                                         
C                                                                               
C     ANALYSIS                                                                  
C     ASSPROP                                                                   
C     BOLOM                                                                     
C     CHKBOLOM                                                                  
C     CHKCONV                                                                   
C     CHKFLUX                                                                   
C     CONVERG1                                                                  
C     CONVERG2                                                                  
C     CONVOLVE                                                                  
C     CONV2D                                                                    
C     ETA                                                                       
C     ETAFUN                                                                    
C     FINDINT                                                                   
C     GETBOUT                                                                   
C     GETOMEGA                                                                  
C     GETOPTPR                                                                  
C     GETPROP                                                                   
C     GETSIZES                                                                  
C     GETTAU                                                                    
C     GETETAZP                                                                  
C     IMAGFN                                                                    
C     MIE                                                                       
C     PSFN                                                                      
C     SETUPETA                                                                  
C     SIZEDIST                                                                  
C     SPECTRAL                                                                  
C     SPFEATUR                                                                  
C     PHILAM                                                                    
C     VISIBILI                                                                  
C     VISI2D                                                                    
c     =======================================================================       
                                                                                

      SUBROUTINE ANALYSIS(model,ETAzp,error)                                    
c     =======================================================================       
c     This subroutine analyzes the solution. It finds the flux conservation         
c     accuracy and evaluates many output quantites like QF(y), TAUF(y),Psi, F1,     
c     the rad.pressure force, dynamical quantities etc.                             
c                                                       [ZI,Mar'96;MN,Mar'99]       
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      DOUBLE PRECISION ugas(npY), qF(npY), vrat(npG,npY), Gamma(npY),           
     &       I1, I2, I3, CMdot, Cve, CM, Cr1                                    
      COMMON /dyn/ ugas, qF, vrat, Gamma, I1, I2, I3, CMdot, Cve, CM,           
     &       Cr1                                                                
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER iL, iY, model, iP, error                                          
      DOUBLE PRECISION ETAzp(npP,npY), QpTd(npG,npY), QpStar(npY), Psi,         
     &        qaux(npL), qaux2(npL), resaux, xP, Planck, QUtot1, r_gd,          
     &        Eps1, aux, deltauF, C1,C2,C3, theta1, ugas_out, s4,Pi,mx          

      Pi = 2.0*ASIN(1.0)                                                        
      r_gd = 200.0                                                              
c     make sure that grids correspond to accepted solution                      
      nY = nYok                                                                 
      DO iY = 1, nY                                                             
        Y(iY) = Yok(iY)                                                         
      END DO                                                                    
      nP = nPok                                                                 
      DO iP = 1, nP                                                             
        P(iP) = Pok(iP)                                                         
      END DO                                                                    

c     spectrum (flux at the outer edge as a function of wavelength)             
      DO iL = 1, nL                                                             
        Spectrum(iL) = dabs(ftot(iL,nY))                                        
c       added in version dusty16.for - to prevent taking log from zero          
c       in Spectral [MN]:                                                       
        IF (Spectrum(iL).LE.1.0D-20) Spectrum(iL) = 1.0D-20                     
      END DO                                                                    

c     analyze bolometric flux error (1/2 of the max spread of fbol)             
      CALL FindErr(fbol,maxFerr,nY)                                             

c     find the flux averaged optical depth, tauF(y)                             
      IF (denstyp.NE.0) THEN                                                    
c      for spherical shell                                                      
       tauF(1) = 0.0                                                            
       DO iY = 2, nY                                                            
c        generate auxiliary function for integration:                           
c        loop over iL (wavelength)                                              
         DO iL = 1, nL                                                          
          qaux(iL)=TAUtot(iL)*ETAzp(1,iY)*dabs(ftot(iL,iY))/lambda(iL)          
         END DO                                                                 
         CALL Simpson(npL,1,nL,lambda,qaux,resaux)                              
         aux = resaux / ETAzp(1,iY)                                             
         deltauF = aux * (ETAzp(1,iY) - ETAzp(1,iY-1))                          
         tauF(iY) = tauF(iY-1) + deltauF                                        
       END DO                                                                   
      ELSE                                                                      
c      for slab                                                                 
       tauF(1) = 0.0                                                            
       DO iY = 1, nY                                                            
c        generate auxiliary function for integration:                           
c        loop over iL (wavelength)                                              
         DO iL = 1, nL                                                          
           qaux(iL)=TAUslb(iL,iY)*dabs(fTot(iL,iY))/lambda(iL)                  
           CALL Simpson(npL,1,nL,lambda,qaux,resaux)                            
           tauF(iY) = resaux                                                    
         END DO                                                                 
       END DO                                                                   
      END IF                                                                    

c     ratio of gravitational to radiation pressure force (isotropic             
c     scattering) per unit volume                                               
c     s4 = (L4sol/Msol)/(4*Pi*G*c*rho_s)/1e-6;                                  
c     rho_s=3000 kg.m-3, grain radius 'a' is in microns, aveV=4/3*Pi*<a^3>      
      IF(denstyp.NE.0) THEN                                                     
      s4 = 1.925 / (4.0*Pi*6.67d-11*3.0d08*3000.0*1.0D-06)                      
c      in case of sigma's from a file aveV=1 (initialized in GetOptPr)
       DO iY = 1, nY                                                            
        DO iL = 1, nL                                                           
          qaux(iL)=(SigmaA(1,iL)+SigmaS(1,iL))/aveV *                           
     &                             dabs(ftot(iL,iY))/lambda(iL)                 
        END DO                                                                  
        CALL Simpson(npL,1,nL,lambda,qaux,resaux)                               
        rg(1,iY) = s4 * resaux / r_gd                                           
c       If dust drift (dynamics case):                                          
        IF (RDW) rg(1,iY) = rg(1,iY)*vrat(1,iY)                                 
       END DO                                                                   
      END IF                                                                    

c     find the Planck averaged absorption efficiencies                          
      DO iY = 1, nY                                                             
c     generate auxiliary function for integration over wavelengths:             
       DO iL = 1, nL                                                            
         qaux(iL) = SigmaA(1,iL) * Us(iL,iY) / lambda(iL)                       
         xP = 14400.0 / Td(1,iY) / lambda(iL)                                   
         qaux2(iL) = SigmaA(1,iL) * Planck(xP) / lambda (iL)                    
       END DO                                                                   
       CALL Simpson(npL,1,nL,lambda,qaux,resaux)                                
       QpStar(iY) = resaux                                                      
       CALL Simpson(npL,1,nL,lambda,qaux2,resaux)                               
       QpTd(1,iY) = resaux                                                      
      END DO                                                                    

c     find parameter Psi (see Ivezic & Elitzur, 1996)                           
c     generate auxiliary function for integration:                              
c     loop over iL (wavelength)                                                 
      DO iL = 1, nL                                                             
        qaux(iL) = SigmaA(1,iL) * Utot(iL,1) / lambda (iL)                      
      END DO                                                                    
      CALL Simpson(npL,1,nL,lambda,qaux,resaux)                                 
      QUtot1 = resaux                                                           
      Psi = QUtot1 / QpTd(1,1)                                                  

      IF(denstyp.NE.0) THEN                                                     
c      ratio r1/r* (see Ivezic & Elitzur, 1996, eq. 27)                         
       r1rs = 0.5 * dsqrt(Psi) * (Tstar / Tsub(1))**2.0                         
      END IF                                                                    

c     Find epsilon - the relative contribution of the diffuse radiation         
      DO iY = 1, nY                                                             
       aux = QpStar(iY)/QpTd(1,iY)/Psi*(Tsub(1)/Td(1,iY))**4.                   
       IF (denstyp.NE.0) aux = aux/ Y(iY)/Y(iY)                                 
       Eps(iY) = 1. - aux                                                       
      END DO                                                                    
      Eps1 = 1.0 - QpStar(1) / QUtot1                                           
c     store these parameters in the storage array                               
      SmC(1,model) = Psi                                                        
      SmC(2,model) = Eps1                                                       
      SmC(3,model) = QpStar(1)                                                  
      SmC(4,model) = QpTd(1,1)                                                  
      SmC(5,model) = maxFerr                                                    

c     additional output quantities                                              
c     bolometric flux at r1 (in W/m2)                                           
c     The constant is 4*Sigma*T^4 (2.27E5 = 4*5.67E-08*(1000**4))               
      F1 = 2.27E5 / Psi * (Tsub(1)/1000.)**4.0                                  
c     inner radius (in cm)                                                      
      Cr1 = 5.53E16 / dsqrt(F1)                                                 
      IF (denstyp.NE.0) THEN                                                    
c      angular diameter of inner cavity if Fbol=1E-6 W/m2                       
       theta1 = 412.6 / dsqrt(F1)                                               
c      check if the pt.source assumption is still obeyed                        
c      (only for BB-type spectrum including EM-function)
       IF(startyp(1).eq.1.OR.startyp(1).eq.2) THEN                              
         mx = sqrt(sqrt(F1/5.67E-08))                                           
         Te_min = 2. * DMAX1(Tsub(1), mx)                                       
       END IF                                                                   
      END IF                                                                    
      IF (denstyp.EQ.0) THEN                                                    
c       Teff for the left illuminating source in slab geometry                  
c       Teff = (F1/sigma)^0.25                                                  
        SmC(7,model) = sqrt(sqrt(F1/5.67E-08))                                  
        IF (ksi.GT.0.) THEN                                                     
c       Teff for the right illuminating source in slab geometry                 
           SmC(8,model) = SmC(7,model)*sqrt(sqrt(ksi))                          
        ELSE                                                                    
           SmC(8,model) = 0.                                                    
        END IF                                                                  
      END IF                                                                    
c     calculate conversion constants for dynamics                               
      IF (denstyp.eq.4.OR.RDW) THEN                                    
        IF (denstyp.EQ.4) THEN                                                  
          I1 = 2. * (1.-pow)/(1.+pow)/tauF(nY)                                  
          I2 = I1                                                               
          I3 = I1 * tauF(nY) / TAUfid                                           
          Gamma(nY) = 0.5                                                       
        END IF                                                                  
c       terminal expansion velocity, full formula:                              
        ugas_out = tauF(nY) * (1.-Gamma(nY)) / (1.-pow)                         
c       The coefficients come from the units conversion                         
        C1 = 0.2845*TAUfid*sqrt(Psi)/I2/(SigExfid/aveV)*                        
     &                                        1.0D06/Tsub(1)/Tsub(1)            
        C2 = 2.040*ugas_out                                                     
        C3 = 6.628*I3*SigExfid/aveV*Gamma(nY)/I1                                
c       from version 2.0 stellar mass is defined as the maximal stellar         
c       mass which does not quench the wind; the calculation is done            
c       with half that mass since any smaller mass will have no effect          
c       on the radial velocity and density profile (see IE '99)                 
c       n.b. Gamma(nY) is removed                                               
        C3 = 6.628*I3*SigExfid/aveV/I1                                          
c       new definitions for output                                              
c       mass-loss rate in Msol/yr                                               
        CMdot = 1.0E-05 * sqrt(C1)                                              
c       terminal expansion velocity in km/s                                     
        Cve = 10.* C2 / sqrt(C1)                                                
c       stellar mass                                                            
        CM = C3                                                                 
      END IF                                                                    

      RETURN                                                                    
      END                                                                       
                                                                                

      SUBROUTINE AssProp(NN,nin,kin,nout,kout)                                  
c     =======================================================================       
c     This subroutine copies arrays nin and kin to nout and kout,                   
c     respectively.                                        [Z.I., Nov. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER NN, i                                                             
      DOUBLE PRECISION nin(NN), kin(NN), nout(NN), kout(NN)                     

      DO i = 1, NN                                                              
        nout(i) = nin(i)                                                        
        kout(i) = kin(i)                                                        
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE BOLOM(q,qBol)                                                  
c     =======================================================================       
c     This subroutine integrates given radiation field, q (function of              
c     wavelength and radial coordinate), over wavelength. q is a matrix             
c     of physical size (npL,npY) [coming from paramet.inc] and real size            
c     (nL,nY) [coming from grids.inc], and qBol is an array of physical size        
c     (npY) and real size nY.                              [Z.I., Mar. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iL, iY                                                            
      DOUBLE PRECISION q(npL,npY), qaux(npL), qBol(npY), resaux                 

c    loop over iY (radial coordinate)                                           
      DO iY = 1, nY                                                             
c       generate auxiliary function for integration                             
c       loop over iL (wavelength)                                               
        DO iL = 1, nL                                                           
          qaux(iL) = q(iL,iY) / lambda (iL)                                     
        END DO                                                                  
        CALL Simpson(npL,1,nL,lambda,qaux,resaux)                               
        qBol(iY) = resaux                                                       
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE ChkBolom(qbol,accur,dev,FbolOK)                                
c     =======================================================================       
c     This subroutine checks if any element of qbol(i), i=1,nY differs for          
c     more than accuracy from 1.0. If so FbolOK = 0, otherwise FbolOK = 1.          
c     dev is maximal deviation from 1.0.                   [Z.I., Mar. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iY, FbolOK                                                        
      DOUBLE PRECISION qBol(npY), accur, dev                                    

      FbolOK = 1                                                                
      dev = 0.0                                                                 
c     loop over iY (radial coordinate)                                          
      DO iY = 1, nY                                                             
        IF (abs(1.0-qBol(iY)).GT.accur) FbolOK = 0                              
        IF (abs(1.0-qBol(iY)).GT.dev) dev = abs(1.0-qBol(iY))                   
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE ChkConv(nY,accuracy,Aold,Anew,Aconv)                           
c     =======================================================================       
c     This subroutine checks convergence of an array A(nY) between values           
c     given in Aold and Anew. If the relative difference for EVERY element          
c     is smaller than accuracy, Aconv is assigned 1, otherwise 0.                   
c                                                          [Z.I., Jul. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, iY, Aconv                                                     
      DOUBLE PRECISION accuracy, Aold(npY), Anew(npY), delta                    

      Aconv = 1                                                                 
c     loop over radial positions                                                
      DO iY = 1, nY                                                             
c       find relative difference                                                
        delta = dabs(Anew(iY)-Aold(iY))                                         
        IF (delta.GT.dabs(Anew(iY))*accuracy) Aconv = 0                         
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE ChkFlux(flux,tolern,consfl,error,ETAzp)                        
c     =======================================================================       
c     Checks the bolometric flux conservation at any point of a given Ygrid.        
c     In case of nonconservation increases the number of points at certain          
c     places. The current criterion is increasing the flux difference from          
c     tolern to its maximum value.                         [MN & ZI,July'96]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER iYins(npY), k, kins, iY, consfl, flag, error, istop, iDm          
      DOUBLE PRECISION Yins(npY),flux(npY), tolern, ff, Yloc,                   
     &       delTAUMax, devfac, devmax, ffold, EtaTemp(npY),ee,                 
     &       ETA, ETAzp(npP,npY)                                                

c     save old grid and values of Eta (important for denstyp = 5 or 6)          
      IF (RDW) THEN                                                             
        DO iY = 1, nY                                                           
          Yprev(iY) = Y(iY)                                                     
          EtaTemp(iY) = ETAdiscr(iY)                                            
        END DO                                                                  
        nYprev = nY                                                             
      END IF                                                                    
      error = 0                                                                 
      kins = 0                                                                  
      devmax = 0.0                                                              
c     maximal delTAU is no more than 2 times the average value                  
      delTAUmax = 2.0*TAUtot(1)*ETAzp(1,nY)/nY                                  
c     maximal deviation from 1.0                                                
      DO iY = 2, nY                                                             
        IF (abs(flux(iY)-1.0).GT.devmax) devmax = abs(flux(iY)-1.0)             
      END DO                                                                    
      ff = 0.0                                                                  
      istop = 0                                                                 
      devfac = 0.1                                                              
c     search for places to improve the grid                                     
      DO WHILE (istop.NE.1)                                                     
        DO iY = 2, nY                                                           
          ffold = ff                                                            
          ff = abs(flux(iY)-1.0)                                                
          flag = 0                                                              
c         if any of these criteria is satisfied insert a point:                 
c         1) if error is increasing too fast                                    
          IF (abs(ff-ffold).GT.devfac*devmax) flag = 1                          
c         2) if delTAU is too large                                             
          IF (TAUtot(1)*(ETAzp(1,iY)-ETAzp(1,iY-1)).GT.                         
     &                             delTAUmax) flag = 1                          
          IF(flag.EQ.1.AND.devmax.GE.tolern) THEN                               
            kins = kins + 1                                                     
            Yins(kins) = Y(iY-1)+0.5*(Y(iY)-Y(iY-1))                            
            iYins(kins) = iY-1                                                  
          END IF                                                                
        END DO                                                                  
        IF (devmax.LT.tolern.OR.devfac.LT.0.01) THEN                            
          istop = 1                                                             
          ELSE                                                                  
          IF (kins.GT.0) istop = 1                                              
        END IF                                                                  
        devfac = devfac / 2.0                                                   
      END DO                                                                    
      IF (kins.EQ.0) THEN                                                       
        IF (consfl.NE.5) consfl = 1                                             
        ELSE                                                                    
c       Add all new points to Y(nY). This gives the new Y(nY+kins).             
c       However, check if npY is large enough to insert all points:             
        IF ((nY+kins).GT.npY) THEN                                              
c        consfl.EQ.5 is a signal that Chkflux was called from SetGrids,         
c        in this case continue without inserting new points. If this is         
c        full problem then give it up.                                          
         IF (consfl.NE.5) THEN                                                  
           consfl = 1                                                           
           ELSE                                                                 
           consfl = 7                                                           
           goto 777                                                             
         END IF                                                                 
         IF (iX.GE.1) THEN                                                      
         write(18,*)' ****************     WARNING   ******************'        
         write(18,*)'  The new Y grid can not accomodate more points!'          
         write(18,'(a,i3)')'   Specified accuracy would require',nY+kins        
         write(18,'(a,i3,a)')'   points, while npY =',npY,'.'                   
         write(18,*)'  For the required accuracy npY must be increased,'        
         write(18,*)'  (see the manual S3.5 Numerical Accuracy).'               
         write(18,*)' *************************************************'        
         END IF                                                                 
         kins = npY - nY                                                        
         iWARNING = iWARNING + 1                                                
         error = 2                                                              
        END IF                                                                  
        DO k = 1, kins                                                          
          CALL SHIFT(Y,npY,nY+k-1,Yins(k),iYins(k)+k-1)                         
        END DO                                                                  
      END IF                                                                    
c     new size of the Y grid                                                    
      nY = nY + kins                                                            
c     intepolate ETAdiscr to new Y grid for denstyp = 5 or 6                    
      DO iY = 1, nY                                                             
        Yloc = Y(iY)                                                            
        IF (iterETA.GT.1) THEN                                                  
          CALL LinInter(npY,nYprev,Yprev,EtaTemp,Yloc,iDm,ee)                   
          ETAdiscr(iY) = ee                                                     
          ELSE                                                                  
          ETAdiscr(iY) = ETA(Yloc)                                              
        END IF                                                                  
      END DO                                                                    

777   RETURN                                                                    
      END                                                                       


      SUBROUTINE Converg1(nY,accuracy,dynrange,Aold,Anew,Aconv,dmax)            
c     =======================================================================       
c     This subroutine checks convergence of an array A(nL,nY) between values        
c     given in Aold and Anew, when the values are larger than dynrange. If          
c     the maximum relative difference is smaller than the required accuracy,        
c     Aconv is assigned 1, otherwise 0.              [Z.I.Jul 96;M.N.Apr.97]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, iY, Aconv                                                     
      DOUBLE PRECISION accuracy, dynrange, Aold(npY), Anew(npY), delta,         
     &       dmax                                                               

      Aconv = 1                                                                 
      dmax = 0.0                                                                
c     loop over radial positions                                                
      DO iY = 1, nY                                                             
c       do it only for elements larger than dynrange                            
        IF (Anew(iY).GE.dynrange) THEN                                          
c           find relative difference                                            
            delta = dabs((Anew(iY)-Aold(iY))/Anew(iY))                          
            IF (delta.GT.dmax) dmax = delta                                     
        END IF                                                                  
      END DO                                                                    
      IF (dmax.GT.accuracy) Aconv = 0                                           

      RETURN                                                                    
      END                                                                       


      SUBROUTINE Converg2(nY,nL,accuracy,dynrange,Aold,Anew,Aconv,dmax)         
c     =======================================================================       
c     This subroutine checks convergence of an array A(nL,nY) between values        
c     given in Aold and Anew, when the values are larger than dynrange. If          
c     the maximum relative difference is smaller than required accuracy,            
c     Aconv is assigned 1, otherwise 0.             [Z.I.Jul 96; M.N.Apr.97]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nL, iY, iL, Aconv                                             
      DOUBLE PRECISION accuracy, dynrange, Aold(npL,npY), Anew(npL,npY),        
     &       delta, dmax                                                        

      Aconv = 1                                                                 
      dmax = 0.0                                                                
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c       loop over radial positions                                              
        DO iY = 1, nY                                                           
c         do it only for elements larger than dynrange                          
          IF (Anew(iL,iY).GE.dynrange) THEN                                     
c           find relative difference                                            
            delta = dabs((Anew(iL,iY)-Aold(iL,iY))/Anew(iL,iY))                 
            IF (delta.GT.dmax) dmax = delta                                     
          END IF                                                                
        END DO                                                                  
      END DO                                                                    
      IF (dmax.GT.accuracy) Aconv = 0                                           

      RETURN                                                                    
      END                                                                       


      SUBROUTINE Convolve(IntOut)                                               
c     =======================================================================       
c     This subroutine convolves intensity IntOut with the point spread              
c     function to produce convolved images ConvInt. The work horse is               
c     subroutine Conv2D, and this subroutine is used to prepare everything.         
c                                                          [Z.I., Jan. 1997]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         
      INTEGER i, j                                                              
      DOUBLE PRECISION IntOut(20,npP+2), yang(npP+2), Youtang, deltaOff,        
     &       Int1D(npP+2),  ConvS, Conv(npP+2), FWHM1max, FWHM2max, PSFN        

c     find the largest FWHMs                                                    
      FWHM1max = FWHM1(1)                                                       
      FWHM2max = FWHM2(1)                                                       
      IF (psftype.LT.3) THEN                                                    
        DO i = 1, NlambdaOut                                                    
          IF (FWHM1(i).GT.FWHM1max) FWHM1max = FWHM1(i)                         
          IF (FWHM2(i).GT.FWHM2max) FWHM2max = FWHM2(i)                         
        END DO                                                                  
      END IF                                                                    
c     scale angular coordinate to theta1                                        
      DO i = 1, nP+2                                                            
        yang(i) = bOut(i) * Theta1 / 2.                                         
      END DO                                                                    
c     generate off-set grid                                                     
      Youtang = Y(nY) * Theta1                                                  
      IF (Youtang.GT.FWHM1max.AND.Youtang.GT.FWHM2max) THEN                     
c       the envelope is well resolved, take impact parameter grid               
        Nconv = nP + 2                                                          
        DO i = 1, Nconv                                                         
          Offset(i) = yang(i)                                                   
        END DO                                                                  
      ELSE                                                                      
c       the envelope is not well resolved, take equidistant grid                
c       to 2FWHM1max, i.e. image will be more or less the PSF itself            
        Nconv = 30                                                              
        deltaOff = 2.0 * FWHM1max / (Nconv-1)                                   
        IF (FWHM2max.GT.FWHM1max) THEN                                          
          deltaOff = 2.0 *FWHM2max / (Nconv-1)                                  
        END IF                                                                  
        DO i = 1, Nconv                                                         
          Offset(i) = deltaOff * 1.0*(i-1)                                      
        END DO                                                                  
      END IF                                                                    
c     convolve intensity wavelength by wavelength                               
      DO j = 1, NlambdaOut                                                      
c       needed to specify wavelength in PSFN                                    
        iLambda = j                                                             
c       generate 1D intensity vector for subrutine Conv2D                       
c       take only diffuse emission, stellar contribution will                   
c       be added below (a shortcut to avoid inaccuracies or too many            
c       points in Conv2D)                                                       
        DO i = 1, nP+2                                                          
          IF (i.LE.2) THEN                                                      
            Int1D(i) = IntOut(j,3)                                              
          ELSE                                                                  
            Int1D(i) = IntOut(j,i)                                              
          END IF                                                                
          CALL CHKRANGE(dynrange,Int1D(i))                                      
          IF (Int1D(i).LT.dynrange) Int1D(i)=0.0                                
        END DO                                                                  
c       convolve                                                                
        CALL Conv2D(npP+2,nP+2,yang,Int1D,1000,NConv,Offset,Conv)               
c       add stellar contribution                                                
        DO i = 1, nP+2                                                          
          ConvS=2.*ASIN(1.0)*(yang(2)**2.)*IntOut(j,1)*PSFN(Offset(i))          
          Conv(i) = Conv(i) + ConvS                                             
        END DO                                                                  
c       scale to 1 at the center                                                
        CALL ScaleTo1(1000,Nconv,Conv)                                          
c       copy 1D convolved intensity to ConvInt                                  
        DO i = 1, Nconv                                                         
          CALL CHKRANGE(dynrange,Conv(i))                                       
          ConvInt(j,i) = Conv(i)                                                
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE Conv2D(NinMax,Nin,Xin,Yin,Noutmax,Nout,Xout,Yout)              
c     =======================================================================       
c     This subroutine convolves intensity Yin(Xin[i]), i=1,Nin with                 
c     the point spread function PSFN(x) (provided as a separate function).          
c     It is assumed that both the intensity yin and PSFN(x) are circularly          
c     symmetric functions of radial coordinate x, i.e., this subroutine             
c     performs two-dimensional convolution. Convolved intensity, Yout, is           
c     evaluated for very position Xout[i], i=1,Nout, as:                            
c            Yout(Xout) = Int[Yin(yloc)*PSF(xloc)*yloc*dyloc*dphi]                  
c     where xloc = sqrt(yloc **2+Xout**2-2*yloc*Xout*cos(phi), with yloc            
c     and phi being dummy integration variables. Declared size of Xin is            
c     NinMax, the one for Xout is NoutMax. The radial integration is done           
c     using subroutine ROMBY and angular integration is done by using               
c     Simpson rule.                                        [Z.I., Jan. 1997]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER NinMax, Nin, NoutMax, Nout, iphi, iXin, Nphi, iXOut               
      DOUBLE PRECISION Xin(NinMax), Yin(NinMax), Xout(NoutMax), A, B,           
     &       Yout(NoutMax), dphi, phi(1000), fphi(1000), int1, int2,            
     &       imagfn                                                             
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         
      INTEGER ftype                                                             
      DOUBLE PRECISION Ckn, Cxout, Cphi, Cqtheta                                
      COMMON /imfn1/ ftype                                                      
      COMMON /imfn2/ Ckn, Cxout, Cphi, Cqtheta                                  
      EXTERNAL imagfn                                                           

c     Parameters for integration:                                               
c     number of angular points                                                  
      Nphi = 9                                                                  
c     step in angle phi                                                         
      dphi = 2.0*ASIN(1.0) / (Nphi-1)                                           
c     flag for imgfn                                                            
      ftype = 1                                                                 
c     Start integrations                                                        
c     loop over output positions                                                
      DO iXout = 1, Nout                                                        
        Cxout = Xout(iXout)                                                     
c       loop over angular wedges (phi integration)                              
        DO iphi = 1, Nphi                                                       
          phi(iphi) = dphi*1.0*(iphi-1)                                         
          Cphi = phi(iphi)                                                      
          fphi(iphi) = 0.0                                                      
c         loop over input radial positions (radial integration)                 
          DO iXin = 1, Nin-1                                                    
            Ckn = 1.0                                                           
            CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int1)                       
            Ckn = 2.0                                                           
            CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int2)                       
c           contribution from this annulus (lin. approx. for intensity)         
            A = Xin(iXin+1)*Yin(iXin) - Xin(iXin)*Yin(iXin+1)                   
            A = A / (Xin(iXin+1)-Xin(iXin))                                     
            B = (Yin(iXin+1)-Yin(iXin)) / (Xin(iXin+1)-Xin(iXin))               
            fphi(iphi) = fphi(iphi) + A*int1 + B*int2                           
          END DO                                                                
        END DO                                                                  
        CALL Simpson(1000,1,Nphi,phi,fphi,Yout(iXout))                          
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      DOUBLE PRECISION FUNCTION ETA(Yy)                                         
c     =======================================================================       
c     This subroutine evaluates the normalized density profile. denstype is         
c     a type of density law: 1 and 2 - power law, 3 for exponential density law,    
c     4 for radiatively driven winds (the gray-body approximation), 5,6 - RDW and   
c     7 - d.d.from a file. pow is parameter describing the choosen density law:     
c     power for 1 and 2, v1/v8 for RDW (the ratio of expansion velocities at the    
c     inner and outer radii), sigma for 3 [i.e. rho = dexp(-(y/sigma)**2)].         
c     Yout is the relative thickness, Yout=rout/r1. Y is the radial position.       
c                                                              [ZI'95; ZI'99]       
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER i, istop, iYdummy                                                 
      DOUBLE PRECISION Yy, C, IntAux, powaux, Prod, eps, factEta                

c     this is an adjustable value regulating the initial Eta approximation for d
      factEta = 0.5                                                             
      IF (Yy.GT.Yout) Yy = Yout                                                 
      IF (iterETA.GT.1) THEN                                                    
c     if this is iteration over ETA (but not the first one) in the case         
c     of dynamical calculation for radiatively driven winds calculate           
c     ETA by linear interpolation of ETA from the previous iteration            
      CALL LinInter(npY,nYprev,Yprev,ETAdiscr,Yy,iYdummy,ETA)                   
c     otherwise use prescribed formulae for different cases                     
      ELSE                                                                      
c     smooth power-law                                                          
      IF (denstyp.EQ.1) THEN                                                    
c       find normalization constant                                             
        IF (pow.NE.1.0) THEN                                                    
          C = (1.0 - Yout**(1.0-pow)) / (pow - 1.0)                             
          ELSE                                                                  
          C = dlog(Yout)                                                        
        ENDIF                                                                   
        IF (Ntr.GE.1) THEN                                                      
          DO i = 1, Ntr                                                         
            powaux = pow - ptr(i)                                               
            IF (powaux.NE.1.0) THEN                                             
              IntAux = (1.0 - Yout**(1.0-powaux)) / (powaux - 1.0)              
            ELSE                                                                
              IntAux = dlog(Yout)                                               
            END IF                                                              
            C = C + IntAux / Ytr(i)**ptr(i)                                     
          END DO                                                                
        ENDIF                                                                   
        C = 1.0 / C                                                             
c       calculate density                                                       
        IF (Yy.GE.(1.0-1.0E-8)) THEN                                            
          ETA = C / Yy**pow                                                     
          IF (Ntr.GE.1) THEN                                                    
            DO i = 1, Ntr                                                       
              ETA = ETA + C * Yy**(ptr(i)-pow) / Ytr(i)**ptr(i)                 
            END DO                                                              
          END IF                                                                
          ELSE                                                                  
          ETA = 0.0                                                             
        ENDIF                                                                   
      END IF                                                                    
c     broken power-law                                                          
      IF (denstyp.EQ.2) THEN                                                    
        Ytr(Ntr+1) = Yout                                                       
c       find normalization constants                                            
        IF (pow.NE.1.0) THEN                                                    
          C = (1.0 - Ytr(1)**(1.0-pow)) / (pow - 1.0)                           
          ELSE                                                                  
          C = dlog(Ytr(1))                                                      
        ENDIF                                                                   
        IF (Ntr.GE.1) THEN                                                      
          DO i = 1, Ntr                                                         
            CALL DoProduct(10,Ytr,ptr,pow,i,Prod)                                 
            IF (ptr(i).NE.1.0) THEN                                             
              IntAux = Ytr(i)**(1.0-ptr(i)) - Ytr(i+1)**(1.0-ptr(i))            
              IntAux = Prod * IntAux / (ptr(i) - 1.0)                           
              ELSE                                                              
              IntAux = Prod * dlog(Ytr(i+1)/Ytr(i))                             
            END IF                                                              
            C = C + IntAux                                                      
          END DO                                                                
        ENDIF                                                                   
        C = 1.0 / C                                                             
c       calculate density                                                       
        IF (Yy.GE.1.0-1.0E-8) THEN                                              
          IF (Yy.LE.Ytr(1)) THEN                                                
            ETA = C / Yy**pow                                                   
          ELSE                                                                  
            istop = 0                                                           
            i = 0                                                               
            DO WHILE (istop.NE.1)                                               
              i = i + 1                                                         
              IF (Yy.LE.Ytr(i+1)) istop = 1                                     
            END DO                                                              
            CALL DoProduct(10,Ytr,ptr,pow,i,Prod)                                 
            ETA = C * Prod / Yy**ptr(i)                                         
          END IF                                                                
        ELSE                                                                    
          ETA = 0.0                                                             
        ENDIF                                                                   
      END IF                                                                    
c     exponential law                                                           
      IF (denstyp.EQ.3) THEN                                                    
        IF (Yy.GE.1.0-1.0E-8) THEN                                              
          ETA = (Yout-1.) * (1.-exp(-pow)) / pow                                
          ETA = exp(-pow*(Yy-1.)/(Yout-1.)) / ETA                               
        ELSE                                                                    
          ETA = 0.0                                                             
        END IF                                                                  
      END IF                                                                    
c     radiatively driven winds (the gray-body approximation)                    
      IF (denstyp.EQ.4) THEN                                                    
        eps = pow                                                               
        IF (Yy.GE.1.0-1.0E-8) THEN                                              
          ETA = (1.+eps)/2./Yy/Yy/sqrt(1.-(1.-eps*eps)/Yy)                      
        ELSE                                                                    
          ETA = 0.0                                                             
        ENDIF                                                                   
      END IF                                                                    
c     radiatively driven winds (full calculation)                               
       IF (RDW) THEN                                                            
c       if this is the first iteration use analytic approximation               
        IF (iterETA.LT.2) THEN                                                  
c         for 5 eps is pow, for 6 assume eps=0.1                                
          IF (denstyp.EQ.5) THEN                                                
            eps = pow                                                           
          ELSE                                                                  
            eps = 0.1                                                           
          END IF                                                                
          IF (Yy.GE.(1.0-1.0E-8)) THEN                                          
            ETA = (1.+eps)/2./Yy/Yy/sqrt(1.-(1.-eps*eps)/Yy)                    
c           empirical improvement for the initial approximation                 
c           good only for large optical depths, but small ones                  
c           are fast anyway (ZI, May99)                                         
            IF (Yy.LE.2) THEN                                                   
               ETA = ETA / (1. + factEta / Yy**10.)                             
            END IF                                                              
          ELSE                                                                  
            ETA = 0.0                                                           
          ENDIF                                                                 
c       or interpolate from the previous solution                               
        ELSE                                                                    
          CALL LinInter(npY,nYprev,Yprev,ETAdiscr,Yy,iYdummy,ETA)               
        END IF                                                                  
      END IF                                                                    
c     user specified function                                                   
      IF (denstyp.EQ.7) THEN                                                    
        IF (Yy.LT.yEta7(nYEta7)) THEN                                           
          CALL LinInter(npY,nYEta7,yEta7,Eta7,Yy,iYdummy,ETA)                   
        ELSE                                                                    
          ETA = Eta7(nYEta7)                                                    
        END IF                                                                  
      END IF                                                                    
c     done                                                                      
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      DOUBLE PRECISION FUNCTION ETAfun(iL,Y)                                    
c     =======================================================================       
c     This subroutine evaluates the normalized density profile as a function        
c     of position and wavelength (in multigrain case). The MAIN purpose of          
c     this function is to provide connection between the overall, prescribed        
c     density distribution, ETA, and wavelength depedent function ETAfun.           
c     It is NOT FINISHED, this is a trivial case, works only for single             
c     grains.                                              [Z.I., Nov. 1995]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER iL                                                                
      DOUBLE PRECISION Y, eta                                                   

c     this is temporary, works only for single component grains                 
c     I need to add subroutines/functions to take care of temperature           
c     dependent and thus variable condensation radii for different              
c     components                                                                
      ETAfun = ETA(Y)                                                           
c     this is to avoid warning when compiling:                                  
      iL = iL                                                                   

      RETURN                                                                    
      END                                                                       


      SUBROUTINE FindInt(nG,alpha,ETAzp)                                        
c     =======================================================================       
c     This subroutine finds the intensity distribution at outer edge and for        
c     user specified wavelengths lamOut. It also evaluates the angular size         
c     of the stellar disk and adds two impact parameters describing the star        
c     to the P grid, thus producing bOut grid. All intensities are indeed           
c     dimensionless quantities lambda*I_lambda/F1 where I_lambda is real            
c     physical quantity defined as usual and F1 is the bolometric flux at           
c     the dust sublimation radius, r1. For conversion to the physical value         
c     lambda*I_lambda, I_lambda from the program has to be multiplied by F1.        
c     F1 can obtained either as:                                                    
c          1) F1 = 4*sigma*Tsub**4/Psi (IE96, eq. 15),                              
c     where Tsub is sublimation temperature and parameter Psi is given in           
c     *.SUM file; or as:                                                            
c          2) F1 = Fbol/alpha1**2 (IE96, eq. 34)                                    
c     where Fbol is the bolometric flux and alpha1 is the angular size of r1        
c     at any particular distance from the envelope (i.e. both F1 and alpha1         
c     correspond to observed quantities). Also note that                            
c         INT(I_lambda(p)*2Pi*P*dP) = f_lambda                                      
c     where I_lambda is the scaled quantity from the program, P is impact           
c     parameter, and f_lambda is the spectral shape F_lambda/Fbol.                  
c                                                          [Z.I., Aug. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER iYfirst, YPequal, Plast                                           
      DIMENSION iYfirst(npP), YPequal(npP), Plast(npY)                          
      COMMON /Yfirst/ iYfirst, YPequal, Plast                                   
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER iL, nG, k, i, iLout, iLstop, iP, iW, nZ, Nzpt, iZ, izloc          
      DOUBLE PRECISION   qaux(npL), Psi,alpha(npG,npY), resaux, QUtot1,         
     &       QpTsub, xP, pst, stelfact, Istell(npL), Planck, w1,                
     &       IntL, IntR, xx , alb, Ids(npL,npP), Ide(npL,npP),w2, lw12,         
     &       ETAzp(npP,npY), numcorr, ETAzpStar, qaux2(npL), z1, z2,            
     &       delz, zloc, wloc, resint, pT, Tz, Idboth, tzp(100), pUtot,         
     &       Semis(100), Sscat(100), IntETA, palb, palf, alfa, exterm,          
     &       Utotloc, Sstem(100), Sstsc(100), Istem(npL,100), Idfront,          
     &       Istsc(npL,100), delTau, factaux, UtotL, UtotR, ep1,                
     &       tauzp1, tauInf                                                     

c     temporary                                                                 
      IF (nG.GT.1.AND.iX.GE.1) THEN                                             
        write(18,*)' FindInt should be fixed, nG>1 !'                           
        stop                                                                    
      END IF                                                                    
c     find impact parameter tangential to the stellar disk                      
c     first find the Planck averaged absorption efficiencies at Y=1             
      DO iL = 1, nL                                                             
        qaux(iL) = SigmaA(1,iL) * Utot(iL,1) / lambda (iL)                      
        xP = 14400.0 / Tsub(1) / lambda(iL)                                     
        qaux2(iL) = SigmaA(1,iL) * Planck(xP) / lambda (iL)                     
      END DO                                                                    
      CALL Simpson(npL,1,nL,lambda,qaux,resaux)                                 
      QUtot1 = resaux                                                           
      CALL Simpson(npL,1,nL,lambda,qaux2,resaux)                                
      QpTsub = resaux                                                           
c     parameter Psi (see Ivezic & Elitzur, 1996, eq. C4)                        
      Psi = QUtot1 / QpTsub                                                     
c     ratio pst = rstar/rsub (see Ivezic & Elitzur, 1996, eq. 27)               
      pst = 2.0 / dsqrt(Psi) * (Tsub(1) / Tstar)**2.0                           
      IF (pst.GE.0.5) THEN                                                      
         IF (iX.GE.1) THEN                                                      
            write(18,*)' FindInt: specified dust temperature at the '           
            write(18,*)' inner radius results in r*/r1 >= 0.5: '                
            write(18,*)'    r*/r1 =', pst                                       
            write(18,*)' This violates some of Dusty`s assumptions'             
            write(18,*)'  ------  Please consult the manual ------'             
            write(18,*)'  ####  r*/r1 changed by hand to 0.5  ####'             
         END IF                                                                 
         pst = 0.5                                                              
      END IF                                                                    
      stelfact = 1.0 / pst / pst / 3.141593                                     
c     generate bOut, i.e. insert two points such that                           
c     bOut(k)=0.999*pst and bOut(k+1)=1.001*pst                                 
      CALL GetbOut(npP,nP,P,pst,bOut,k)                                         
c     correction for numerical errors in tau                                    
      numcorr = 1. / TAUtot(1)                                                  
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c        stellar intensity, Istell (extinction already included)                
         Istell(iL) = fs(iL,nY) * stelfact                                      
c        total optical depth along a line of sight                              
         tauOut(iL) = numcorr*TAUtot(iL)                                        
      END DO                                                                    
c     generate diffuse intensities, Ide (emission) and Ids (scat)               
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
        DO iP = 1, nP                                                           
c         maximal number of points along tangential position, z                 
          nZ = nY + 1 - iYfirst(iP)                                             
c         starting value for local radius                                       
          IF (P(iP).GE.1.0) THEN                                                
            w2 = P(iP)                                                          
          ELSE                                                                  
            w2 = 1.0                                                            
          END IF                                                                
c         initialize intensities                                                
          Ide(iL,iP) = 0.0                                                      
          Ids(iL,iP) = 0.0                                                      
          IF (iP.LE.k+1) THEN                                                   
            Istem(iL,iP) = 0.0                                                  
            Istsc(iL,iP) = 0.0                                                  
          END IF                                                                
c         total optical depth along this impact parameter                       
          ep1 = ETAzp(iP,nY+1-iYfirst(iP))*TAUtot(iL)                           
c         loop over z, i.e. steps over points crossing the y grid               
          DO iZ = 2, nZ                                                         
c           index for the ending local radius                                   
            iW = iYfirst(iP) + iZ - 1                                           
c           local boundary radii                                                
            w1 = w2                                                             
            w2 = Y(iW)                                                          
c           corresponding displacements along a line of sight                   
            z1 = sqrt(abs(w1**2.-P(iP)**2.))                                    
            z2 = sqrt(abs(w2**2.-P(iP)**2.))                                    
c           # of pts. for z integration, should increase with deltaTau          
c           it is messy because INT function which would do the job is          
c           not in F77 standard set                                             
            Nzpt = 5                                                            
            delTau = (ETAzp(iP,iW)-ETAzp(iP,iW-1))*TAUtot(iL)                   
            IF (delTau.GT.1) Nzpt = 10                                          
            IF (delTau.GT.5) Nzpt = 20                                          
            IF (delTau.GT.10) Nzpt = 30                                         
            IF (delTau.GT.20) Nzpt = 40                                         
            IF (delTau.GT.50) Nzpt = 50                                         
            delz = (z2-z1) / (Nzpt-1)                                           
c           powers for power-law interpolations between 2 y pts.                
            lw12 = dlog(Y(iW-1)/Y(iW))                                          
c           for T                                                               
            pT = dlog(Td(1,iW)/Td(1,iW-1)) / lw12                               
c           for albedo                                                          
            IF (omega(iL,iW-1).GT.0.0.AND.omega(iL,iW).GT.0.0) THEN             
              palb = dlog(omega(iL,iW)/omega(iL,iW-1)) / lw12                   
            ELSE                                                                
              palb = 0.0                                                        
            END IF                                                              
c           for Utot                                                            
            UtotL = Utot(iL,iW-1)                                               
            UtotR = Utot(iL,iW)                                                 
            CALL CHKRANGE(dynrange,UtotL)                                       
            CALL CHKRANGE(dynrange,UtotR)                                       
            IF (UtotL.GT.0.0.AND.UtotR.GT.0) THEN                               
              pUtot = dlog(UtotR/UtotL) / lw12                                  
            ELSE                                                                
              pUtot = 0.0                                                       
            END IF                                                              
c           for alpha                                                           
            palf = dlog(alpha(1,iW)/alpha(1,iW-1)) / lw12                       
c           tauzp between z=0 and z=z1                                          
            tauzp1 = ETAzp(iP,iZ-1)*TAUtot(iL)                                  
c           integrate between adjacent grid points                              
            DO izloc = 1, Nzpt                                                  
              zloc = z1 + (izloc-1)*delz                                        
              wloc = sqrt(zloc**2 + P(iP)**2)                                   
c             find local TAUzp(w(z))-TAUzp(w1=w(z1))                            
              tzp(izloc) = IntETA(P(iP),iW-1,w1,wloc)*TAUtot(iL)                
c             find Tz = T(zloc) = T(wloc), this works for single                
c             size grains only; for multigrain case one needs to                
c             get Semis by summation over all Td                                
              Tz = Td(1,iW-1) * (Y(iW-1)/wloc)**pT                              
              xP = 14400/lambda(iL)/Tz                                          
c             power-law interpolation for albedo                                
              alb = omega(iL,iW-1) * (Y(iW-1)/wloc)**palb                       
c             power-law interpolation for Utot                                  
              IF (UtotL.GT.0) THEN                                              
                UtotLoc = UtotL * (Y(iW-1)/wloc)**pUtot                         
              ELSE                                                              
                UtotLoc = 0.0                                                   
              END IF                                                            
              CALL CHKRANGE(dynrange,UtotLoc)                                   
c             power-law interpolation for alpha                                 
              alfa = alpha(1,iW-1) * (Y(iW-1)/wloc)**palf                       
c             source functions (wloc**2 because D uses scaled quant.)           
              factaux = 1 / wloc**2 / (4 * 3.14159)                             
              Semis(izloc) = (1-alb) * alfa * Planck(xP) * factaux              
              Sscat(izloc) = alb * UtotLoc * factaux                            
c             check for the dynamic range                                       
              CALL CHKRANGE(dynrange,Semis(izloc))                              
              CALL CHKRANGE(dynrange,Sscat(izloc))                              
c             optical depth from infinity along the line of sight               
              tauInf = ep1 - tauzp1 - tzp(izloc)                                
c             for a line of sight terminating on the star find                  
c             contribution only from the front part of the envelope             
              IF (iP.LE.k+1) THEN                                               
                 IF (tauInf.LT.50) THEN                                         
                   exterm = dexp(-tauInf)                                       
                 ELSE                                                           
                   exterm = 0.0                                                 
                 END IF                                                         
                 Sstem(izloc) = Semis(izloc) * exterm                           
                 Sstsc(izloc) = Sscat(izloc) * exterm                           
              END IF                                                            
c             otherwise take both the front and back contributions              
              IF (tauInf.LT.50) THEN                                            
                 exterm = dexp(-tauInf)+dexp(-tauInf-ep1)                       
              ELSE                                                              
                 exterm = 0.0                                                   
              END IF                                                            
              Semis(izloc) = Semis(izloc) * exterm                              
              Sscat(izloc) = Sscat(izloc) * exterm                              
c             end of local loop over z                                          
            END DO                                                              
c           integrate and add contribution from this step                       
            CALL SIMPSON(100,1,Nzpt,tzp,Semis,resint)                           
            CALL CHKRANGE(dynrange,resint)                                      
            Ide(iL,iP) = Ide(iL,iP) + resint                                    
            CALL SIMPSON(100,1,Nzpt,tzp,Sscat,resint)                           
            CALL CHKRANGE(dynrange,resint)                                      
            Ids(iL,iP) = Ids(iL,iP) + resint                                    
            IF (iP.LE.k+1) THEN                                                 
              CALL SIMPSON(100,1,Nzpt,tzp,Sstem,resint)                         
              CALL CHKRANGE(dynrange,resint)                                    
              Istem(iL,iP) = Istem(iL,iP) + resint                              
              CALL SIMPSON(100,1,Nzpt,tzp,Sstsc,resint)                         
              CALL CHKRANGE(dynrange,resint)                                    
              Istsc(iL,iP) = Istsc(iL,iP) + resint                              
            END IF                                                              
c         end of loop over z                                                    
          END DO                                                                
c       end of loop over impact parameter, iP                                   
        END DO                                                                  
c     end of loop over wavelengths, iL                                          
      END DO                                                                    
c     add all intensities, Istell, Ide, Ids                                     
      DO iL = 1, nL                                                             
c       interpolate optical depth  at pstar                                     
        IF (iL.EQ.iLfid) THEN                                                   
          ETAzpStar = (ETAzp(k,nY) - ETAzp(k-1,nY))                             
          ETAzpStar = ETAzpStar * (pst-P(k-1)) / (P(k) - P(k-1))                
          ETAzpStar = ETAzp(k-1,nY) + ETAzpStar                                 
        END IF                                                                  
c       find diffuse contribution at pstar (by linear interpolation)            
        Idfront = Istsc(iL,k)+Istem(iL,k)-Istsc(iL,k-1)-Istem(iL,k-1)           
        Idfront = Idfront * (pst-P(k-1)) / (P(k) - P(k-1))                      
        Idfront = Idfront + Istsc(iL,k-1) + Istem(iL,k-1)                       
        Idboth = Ids(iL,k) + Ide(iL,k) - Ids(iL,k-1) - Ide(iL,k-1)              
        Idboth = Idboth * (pst-P(k-1)) / (P(k) - P(k-1))                        
        Idboth = Idboth + Ids(iL,k-1) + Ide(iL,k-1)                             
c       first for p<pstar, all three contributions                              
        DO i = 1, k-1                                                           
          Intens(iL,i) = Istell(iL) + Istsc(iL,i) + Istem(iL,i)                 
          IF (iL.EQ.iLfid)                                                      
     &        tauZout(i) = ETAzp(i,nY)/ETAzp(1,nY)                              
        END DO                                                                  
c       barely on the stellar disk                                              
        Intens(iL,k) = Istell(iL) + Idfront                                     
        tauZout(k) = ETAzpStar/ETAzp(1,nY)                                      
c       barely off the stellar disk                                             
        Intens(iL,k+1) = Idboth                                                 
        tauZout(k+1) = 2. * tauZout(k)                                          
c       all other p>pstar                                                       
        DO i = k, nP                                                            
          Intens(iL,i+2) = Ids(iL,i)+Ide(iL,i)                                  
          IF (iL.EQ.iLfid) THEN                                                 
            nZ = nY + 1 - iYfirst(i)                                            
            tauZout(i+2) = 2. * ETAzp(i,nZ)/ETAzp(1,nY)                         
          END IF                                                                
        END DO                                                                  
      END DO                                                                    
c     check dynamic range                                                       
      DO iL = 1, nL                                                             
        DO i = 1, nP+2                                                          
          CALL CHKRANGE(dynrange,Intens(iL,i))                                  
        END DO                                                                  
      END DO                                                                    
c     now interpolate Intens(lambda) to lamOut                                  
      DO iLout = 1, NlambdaOut                                                  
c       bracket the needed wavelength                                           
        iLstop = 0                                                              
        iL = 0                                                                  
        DO WHILE (iLstop.EQ.0)                                                  
          iL = iL + 1                                                           
          IF (lambda(iL).GT.LambdaOut(iLout)) iLstop = 1                        
          IF (iL.EQ.nL) iLstop = 1                                              
        END DO                                                                  
c       interpolate intensity                                                   
        xx = (LambdaOut(iLout)-lambda(iL-1))/(lambda(iL)-lambda(iL-1))          
        DO i = 1, nP+2                                                          
          IntL = Intens(iL-1,i)                                                 
          IntR = Intens(iL,i)                                                   
          IntOut(iLout,i) = IntL + xx*(IntR - IntL)                             
          CALL CHKRANGE(dynrange,IntOut(iLout,i))                               
        END DO                                                                  
      END DO                                                                    

999   RETURN                                                                    
      END                                                                       


      SUBROUTINE GetbOut(npP,nP,P,pstar,bOut,k)                                 
c     =======================================================================       
c     This subroutine inserts two impact parameters corresponding to pstar,         
c     producing bOut(nP+2) from P(nP). The inserted elements are bOut(k) and        
c     bOut(k+1)                                            [Z.I., Aug. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npP, nP, k, kstop, i                                              
      DOUBLE PRECISION P(npP), bOut(npP+2), pstar                               

      k = 0                                                                     
      kstop = 0                                                                 
      DO WHILE (kstop.NE.1)                                                     
        k = k + 1                                                               
        bOut(k) = P(k)                                                          
        IF (1.001*pstar.LE.P(k).OR.k.EQ.nP) kstop = 1                           
      END DO                                                                    
      IF (0.999*pstar.GT.P(k-1)) THEN                                           
        bOut(k) = 0.999*pstar                                                   
        ELSE                                                                    
        bOut(k) = 0.5*(P(k-1)+1.001*pstar)                                      
      END IF                                                                    
      IF (1.001*pstar.LT.P(k)) THEN                                             
        bOut(k+1) = 1.001*pstar                                                 
        ELSE                                                                    
        bOut(k+1) = 0.5*(P(k)+0.999*pstar)                                      
      END IF                                                                    
      DO i = k, nP                                                              
        bOut(i+2) = P(i)                                                        
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE getOmega(nG)                                                   
c     =======================================================================       
c     This subroutine generates albedo omega(iL,iY) from the abs/sca cross-         
c     sections and the component abundancies. This is temporary (trivial)           
c     version  for single size grains.                     [Z.I., Mar. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER  iG, nG, iL, iY                                                   

c     generate overall albedo through the envelope                              
c      ** this is for future multigrain code **                                 
c      ** for single grains it is trivial **                                    
      DO iY = 1, nY                                                             
c       calculate albedo                                                        
        DO iL = 1, nL                                                           
          omega(iL,iY) = SigmaS(nG,iL) / (SigmaA(nG,iL) + SigmaS(nG,iL))        
        END DO                                                                  
c       calculate relative abundances                                           
        DO iG = 1, nG                                                           
          abund(iG,iY) = 1.0                                                    
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE getOptPr(nG,nameQ,nameNK,er)                                   
c     =====================================================================         
c     This subroutine calculates the absorption and scattering efficiences          
c     Qabs and Qsca in the wavelength range of the code or in case of               
c     user supplied efficiences reads them from a file.                             
c                                                     [ZI Mar96; MN Aug97]          
c     =====================================================================         
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER npLnk                                                             
      PARAMETER (npLnk=98)                                                      
      DOUBLE PRECISION n_sil_ow(npLnk),k_sil_ow(npLnk),n_sil_oc(npLnk),         
     &       k_sil_oc(npLnk), n_sil_dl(npLnk), k_sil_dl(npLnk),                 
     &       n_amc_hn(npLnk), k_amc_hn(npLnk), n_sic_pg(npLnk),                 
     &       k_sic_pg(npLnk), n_gr1_dl(npLnk), k_gr1_dl(npLnk),                 
     &       n_gr2_dl(npLnk), k_gr2_dl(npLnk), lam_nk(npLnk)                    
      COMMON /nkdat/ n_sil_ow, k_sil_ow, n_sil_oc, k_sil_oc,                    
     &               n_sil_dl, k_sil_dl, n_amc_hn, k_amc_hn,                    
     &               n_sic_pg, k_sic_pg, n_gr1_dl, k_gr1_dl,                    
     &               n_gr2_dl, k_gr2_dl, lam_nk                                 
      CHARACTER*235 nameQ(npG), nameNK(10), Fname, dummy*132                    
      INTEGER iG, nG, io1, iL, nLin, iLaux, Nprop, Na, iA, iC, iCuser,          
     &        er, Nmax, npA                                                     
c     Nmax is the number of records in user supplied file with opt.prop.        
c     and npA is the dimension of the array of grain sizes                      
      PARAMETER (Nmax=10000, npA=100)                                           
      DOUBLE PRECISION aa,bb,cc,lambdain(Nmax),Qain(Nmax),Qsin(Nmax),           
     &       n(npL),k(npL), aQabs(npA,npL),aQsca(npA,npL), amax, Pi,            
     &       nsd(npA), a(npA), faux1(npA), faux2(npA), f(npA), int,             
     &       ala(Nmax), SigAbs(npA,npL), SigSca(npA,npL), SizeDist,             
     &       aQa(Nmax), aQs(Nmax),  Cnorm, a3ave,                               
     &       n_int(npL), k_int(npL)                                             
c ----------------------------------------------------------------              
c     this should never change                                                  
      nL = npL                                                                  
      Nprop = 7                                                                 
c ----------------------------------------------------------------              
      Pi = 2.0*ASIN(1.0)                                                        
      er = 0                                                                    
c     first check that the user supplied wavelength grid is                     
c     monotonously increasing                                                   
      IF (top.LT.3) THEN                                                        
c       calculate efficiencies from n and k by Mie theory                       
c       generate the size array                                                 
        IF (szds.GT.2) THEN                                                     
          amax = 5.0*a2                                                         
        ELSE                                                                    
          amax = a2                                                             
        END IF                                                                  
        IF (dabs(a1-a2).LE.1.d-3) THEN                                          
          nA = 1                                                                
        ELSE                                                                    
          nA =50                                                                
        END IF                                                                  
c       Build-up the array of sizes a(nA)                                       
        CALL GETsizes(npA,nA,a1,amax,a)                                         
c       evaluate the normalization constant for the size                        
c       distribution nsd(nA)                                                    
        DO iA = 1, nA                                                           
          nsd(iA) = SizeDist(qsd,a(iA),szds,a2)                                 
        END DO                                                                  
        CALL PowerInt(npA,1,nA,a,nsd,Cnorm)                                     
c       find the average grain volume aveV (needed in dynamics)                 
        IF(dabs(a1-a2).LE.1.d-3) THEN                                           
          aveV = 4./3.*Pi*a1**3                                                 
        ELSE                                                                    
           DO iA = 1, nA                                                        
             faux1(iA)=nsd(iA)*a(iA)**3                                         
           END DO                                                               
           CALL PowerInt(npA,1,nA,a,faux1,a3ave)                                
           aveV = 4./3.*Pi*a3ave/Cnorm                                          
        END IF                                                                  
c      --  LOOP OVER SUPPORTED COMPONENTS --                                    
        DO iC= 1, Nprop                                                         
          f(iC) = xC(iC)                                                        
c         assign optical properties                                             
          IF (iC.EQ.1) CALL AssProp(npLnk,n_sil_ow,k_sil_ow,n,k)                
          IF (iC.EQ.2) CALL AssProp(npLnk,n_sil_oc,k_sil_oc,n,k)                
          IF (iC.EQ.3) CALL AssProp(npLnk,n_sil_dl,k_sil_dl,n,k)                
          IF (iC.EQ.4) CALL AssProp(npLnk,n_gr1_dl,k_gr1_dl,n,k)                
          IF (iC.EQ.5) CALL AssProp(npLnk,n_gr2_dl,k_gr2_dl,n,k)                
          IF (iC.EQ.6) CALL AssProp(npLnk,n_amc_hn,k_amc_hn,n,k)                
          IF (iC.EQ.7) CALL AssProp(npLnk,n_sic_pg,k_sic_pg,n,k)                
c         interpolate from opt. prop. grid to working grid                      
          DO iL = 1, nL                                                         
             CALL LinInter(npLnk,npLnk,lam_nk,n,lambda(iL),iLaux,aa)            
             n_int(iL) = aa                                                     
             CALL LinInter(npLnk,npLnk,lam_nk,k,lambda(iL),iLaux,aa)            
             k_int(iL) = aa                                                     
          END DO                                                                
c         calculate Qabs and Qsca (on working grid, i.e. lambda)                
          CALL MIE(npL,nL,lambda,n_int,k_int,npA,nA,a,1,aQabs,aQsca)            
c         for each lambda integrate Pi*a^2*Qext with n(a)da                     
          DO iL = 1, nL                                                         
            DO iA = 1, nA                                                       
              faux1(iA)=nsd(iA)*aQabs(iA,iL)*Pi*a(iA)**2                        
              faux2(iA)=nsd(iA)*aQsca(iA,iL)*Pi*a(iA)**2                        
            END DO                                                              
            CALL PowerInt(npA,1,nA,a,faux1,int)                                 
            sigAbs(iC,iL) = int/Cnorm                                           
            CALL PowerInt(npA,1,nA,a,faux2,int)                                 
            sigSca(iC,iL) = int/Cnorm                                           
          END DO                                                                
        END DO                                                                  
        IF (top.EQ.2) THEN                                                      
c         --  LOOP OVER USER SUPPLIED COMPONENTS --                             
          DO iCuser = 1, Nfiles                                                 
            iC = Nprop + iCuser                                                 
            f(iC) = xCuser(iCuser)                                              
c           read in optical properties                                          
            Fname = nameNK(iCuser)                                              
            CALL GetProp(npL,lambda,nL,Fname,n,k,er)                            
            IF (er.EQ.3) goto 999                                               
c           calculate Qabs and Qsca                                             
            CALL MIE(npL,nL,lambda,n,k,npA,nA,a,1,aQabs,aQsca)                  
c           for each lambda integrate Pi*a^2*Qext with n(a)da                   
            DO iL = 1, nL                                                       
              DO iA = 1, nA                                                     
                faux1(iA)=nsd(iA)*aQabs(iA,iL)*Pi*a(iA)**2                      
                faux2(iA)=nsd(iA)*aQsca(iA,iL)*Pi*a(iA)**2                      
              END DO                                                            
              CALL PowerInt(npA,1,nA,a,faux1,int)                               
              sigAbs(iC,iL) = int/Cnorm                                         
              CALL PowerInt(npA,1,nA,a,faux2,int)                               
              sigSca(iC,iL) = int/Cnorm                                         
            END DO                                                              
          END DO                                                                
        ELSE                                                                    
          Nfiles = 0                                                            
        END IF                                                                  
c       mix them together (syntetic grain model)                                
        DO iL = 1, nL                                                           
          SigmaA(1,iL) = 0.0                                                    
          SigmaS(1,iL) = 0.0                                                    
          DO iC= 1, Nprop+Nfiles                                                
            SigmaA(1,iL) = SigmaA(1,iL) + f(iC) * sigAbs(iC,iL)                 
            SigmaS(1,iL) = SigmaS(1,iL) + f(iC) * sigSca(iC,iL)                 
          END DO                                                                
        END DO                                                                  
      ELSE                                                                      
c     this is for top.GE.3                                                      
c      initialize aveV for this case
       aveV = 1.
c       read in lambda grid and optical properties                              
        DO iG = 1, nG                                                           
          open(1,ERR=998,file=nameQ(iG),STATUS='OLD')                           
          read(1,'(a)',ERR=998)dummy                                            
          read(1,'(a)',ERR=998)dummy                                            
          read(1,'(a)',ERR=998)dummy                                            
          iL = 0                                                                
          io1 = 0                                                               
          DO WHILE (io1.GE.0)                                                   
            read(1,*,END=900,ERR=998,iostat=io1) aa, bb, cc                     
            IF (io1.GE.0) THEN                                                  
              iL = iL + 1                                                       
              lambdain(iL) = aa                                                 
              Qain(iL) = bb                                                     
              Qsin(iL) = cc                                                     
            END IF                                                              
          END DO                                                                
900       close(1)                                                              
          IF (iL.LT.2) goto 998                                                 
          nLin = iL                                                             
c         if input wavelengths in descending order turn them around             
          IF (lambdain(1).GT.lambdain(2)) THEN                                  
            DO iL = 1, nLin                                                     
              ala(iL) = lambdain(iL)                                            
              aQa(iL) = Qain(iL)                                                
              aQs(iL) = Qsin(iL)                                                
            END DO                                                              
            DO iL = 1, nLin                                                     
              lambdain(iL) = ala(nLin+1-iL)                                     
              Qain(iL) = aQa(nLin+1-iL)                                         
              Qsin(iL) = aQs(nLin+1-iL)                                         
            END DO                                                              
          END IF                                                                
c         interpolate to Dusty's wavelength grid                                
          DO iL = 1, nL                                                         
            CALL LinInter(Nmax,nLin,lambdain,Qain,lambda(iL),iLaux,aa)          
            SigmaA(iG,iL) = aa                                                  
            CALL LinInter(Nmax,nLin,lambdain,Qsin,lambda(iL),iLaux,aa)          
            SigmaS(iG,iL) = aa                                                  
          END DO                                                                
        END DO                                                                  
      END IF                                                                    
      goto 999                                                                  
998   write(12,*)' ***  FATAL ERROR IN DUSTY  ***********'                      
      write(12,*)' File with optical properties:'                               
      write(12,'(a2,a70)')'  ',nameQ(iG)                                        
      write(12,*)' is missing or not properly formatted?!'                      
      write(12,*)' **************************************'                      
      close(12)                                                                 
      er = 3                                                                    

999   RETURN                                                                    
      END                                                                       


      SUBROUTINE GetProp(npL,lambda,nL,Fname,en,ek,error)                       
c     =======================================================================       
c     This subroutine reads optical properties en(i,j), ek(i,j) from file           
c     fname(Nf), with i=Nf, j=1..NLL(Nf), and interpolates them onto                
c     wavelength grid lambda(1..nL)                        [Z.I., Mar. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      CHARACTER*235 Fname                                                       
      CHARACTER*232 line                                                        
      INTEGER i, nL, iloc, iL, npL, io1, error, Nmax                            
c     Nmax is the number of records in the user supplied file                   
      PARAMETER (Nmax=10000)	                                                   
      DOUBLE PRECISION en(npL), ek(npL), lambda(npL), pw(Nmax),                 
     &       pren(Nmax), pimn(Nmax), a(Nmax), b(Nmax), c(Nmax), aa,             
     &       bb, cc                                                             

      error = 0                                                                 
      open(2,ERR=998,file=Fname,STATUS='OLD')                                   
c     read in a header from the input file                                      
      DO i = 1, 7                                                               
        read(2,'(a)',ERR=998)line                                               
      END DO                                                                    
c     read in input data                                                        
      iL = 0                                                                    
      io1 = 0                                                                   
      DO WHILE (io1.GE.0)                                                       
        read(2,*,END=900,ERR=998,iostat=io1) aa, bb, cc                         
        IF (io1.GE.0) THEN                                                      
          iL = iL + 1                                                           
          pw(iL) = aa                                                           
          pren(iL) = bb                                                         
          pimn(iL) = cc                                                         
        END IF                                                                  
      END DO                                                                    
900   close(2)                                                                  
      IF (iL.LT.2) goto 998                                                     
c     if input wavelengths in descending order turn them around                 
      IF (pw(1).GT.pw(2)) THEN                                                  
        DO i = 1, iL                                                            
          a(i) = pw(i)                                                          
          b(i) = pren(i)                                                        
          c(i) = pimn(i)                                                        
        END DO                                                                  
        DO i = 1, iL                                                            
          pw(i) = a(iL+1-i)                                                     
          pren(i) = b(iL+1-i)                                                   
          pimn(i) = c(iL+1-i)                                                   
        END DO                                                                  
      END IF                                                                    
c     interpolate                                                               
      DO i = 1, nL                                                              
        CALL LinInter(Nmax,iL,pw,pren,lambda(i),iloc,en(i))                     
        CALL LinInter(Nmax,iL,pw,pimn,lambda(i),iloc,ek(i))                     
      END DO                                                                    
      goto 999                                                                  
998   write(12,*)' ***  FATAL ERROR IN DUSTY  ***********'                      
      write(12,*)' File with optical properties:'                               
      write(12,'(a2,a70)')'  ',Fname                                            
      write(12,*)' is missing or not properly formatted?!'                      
      write(12,*)' **************************************'                      
      close(12)                                                                 
      error = 3                                                                 

999   RETURN                                                                    
      END                                                                       


      SUBROUTINE GETsizes(NN,N,x1,x2,x)                                         
c     =======================================================================       
c     This subroutine generates an array x(i=1..N) of physical size NN,             
c     with N elements logarithmically spaced between x1 and x2.                     
c                                                  [ZI,Aug'96;MN,Nov'97]            
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER NN, N, i                                                          
      DOUBLE PRECISION x(NN), x1, x2, fac, pw1, pw2                             

      IF (N.GT.1) THEN                                                          
        pw1 = 1.0/(N-1)                                                         
        fac = (x2/x1)**pw1                                                      
        DO i = 1, N                                                             
          pw2 = 1.0*(i-1)                                                       
          x(i) = x1*fac **pw2                                                   
        END DO                                                                  
      ELSE                                                                      
        x(1) = x1                                                               
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE getTau(model,nG,TAU1,TAU2,TAUin,Nrec,GridType,Nmodel)          
c     =======================================================================       
c     This subroutine generates total optical depth TAUtot.                         
c                                                          [Z.I., Mar. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER model, Nmodel, nG, iL, GridType, Nrec                             
      DOUBLE PRECISION faux1(npL), faux2(npL), SigAfid, SigSfid, q,             
     &       TAUin(Nrec), TAU1, TAU2                                            

      IF (nG.GT.1) THEN                                                         
        write(18,*)'Fix getTAU, nG>1 !'                                         
        stop                                                                    
      END IF                                                                    
      TAUmax = 0.0                                                              
      IF(GridType.EQ.3) THEN                                                    
        TAUfid = TAUin(model)                                                   
      ELSE                                                                      
c      calculate TAUfid for given model                                         
       IF (model.EQ.1) THEN                                                     
        TAUfid = TAU1                                                           
       ELSE                                                                     
        IF (model.EQ.Nmodel) THEN                                               
          TAUfid = TAU2                                                         
        ELSE                                                                    
          IF (GridType.EQ.1) THEN                                               
            q =  (TAU2 - TAU1)/(Nmodel-1.0)                                     
            TAUfid = TAU1 + q*(model-1.0)                                       
          ELSE                                                                  
            q = dexp(log(TAU2/TAU1)/(Nmodel-1.0))                               
            TAUfid = TAU1 * q**(model-1.0)                                      
          END IF                                                                
        END IF                                                                  
       END IF                                                                   
      END IF                                                                    
c     generate TAUtot and find TAUmax                                           
      DO iL = 1, nL                                                             
        faux1(iL) = SigmaA(nG,iL)                                               
        faux2(iL) = SigmaS(nG,iL)                                               
      END DO                                                                    
      IF (lamfid.LT.lambda(1)) THEN                                             
        write(12,*)' Fiducial wavelength was too small.'                        
        write(12,'(a8,e9.3,a17)')' Using ',lambda(1),' micron instead.'         
      END IF                                                                    
      IF (lamfid.GT.lambda(nL)) THEN                                            
        write(12,*)' Fiducial wavelength was too large.'                        
        write(12,'(a8,e9.3,a17)')' Using ',lambda(nL),' micron instead.'        
      END IF                                                                    
      CALL LinInter(npL,nL,lambda,faux1,lamfid,iLfid,SigAfid)                   
      CALL LinInter(npL,nL,lambda,faux2,lamfid,iLfid,SigSfid)                   
c     extinction efficiency at fiducial wavelength                              
      SigExfid = SigAfid + SigSfid                                              
      DO iL = 1, nL                                                             
        TAUtot(iL) = TAUfid*(SigmaA(nG,iL) + SigmaS(nG,iL)) / SigExfid          
        IF (TAUtot(iL).GE.TAUmax) TAUmax = TAUtot(iL)                           
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE getETAzp(ETAzp)                                                
c     =======================================================================       
c     This function calculates ETAzp(iP,iZ) along the line of sight with            
c     impact parameter P(iP) and iZ=1, nZ. Here iZ = 1 corresponds to z=0           
c     and iZ=nZ to the outer edge. Other grid points coincide with the              
c     radial grid. The method used is spline approximation for normalized           
c     density distribution ETA, with subsequent z-integration performed             
c     analytically in function IntETA                                               
c                                                   [ZI,Feb'95; MN,Aug'97]          
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iYfirst, YPequal, Plast                                           
      DIMENSION iYfirst(npP), YPequal(npP), Plast(npY)                          
      COMMON /Yfirst/ iYfirst, YPequal, Plast                                   
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER iP, nZ, iZ, iW                                                    
      DOUBLE PRECISION ETAzp(npP,npY), IntETA, auxEta, w1, w2                   

c     loop over impact parameters                                               
      DO iP = 1, nP                                                             
c       maximal number of points along tangential position, z                   
        nZ = nY + 1 - iYfirst(iP)                                               
c       starting values for z and ETAzp(iP,iZ)                                  
        IF (P(iP).GE.1.0) THEN                                                  
          w2 = P(iP)                                                            
          ELSE                                                                  
          w2 = 1.0                                                              
        END IF                                                                  
c       initialize ETAzp(iP,iZ)*TAUtot(iL)                                      
        ETAzp(iP,1) = 0.0                                                       
c       loop over z                                                             
        DO iZ = 2, nZ                                                           
c         index for local radius, w2                                            
          iW = iYfirst(iP) + iZ - 1                                             
c         limits for integration                                                
          w1 = w2                                                               
          w2 = Y(iW)                                                            
c           find next step in ETAzp                                             
            auxEta = IntETA(P(iP),iW-1,w1,w2)                                   
c           add next step in ETAzp                                              
            ETAzp(iP,iZ) = ETAzp(iP,iZ-1) + auxEta                              
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       
                                                                                

      DOUBLE PRECISION FUNCTION IMAGFN(y)                                       
c     =======================================================================       
c     This function evaluates auxiliary functions needed to produce                 
c     visibility curves and convolved images. It is called from the image           
c     integration subroutine ROMBY.                        [Z.I., Jan. 1997]        
c     =======================================================================       
      IMPLICIT none                                                             
      DOUBLE PRECISION x, y, PSFN, Bessel                                       
      INTEGER ftype                                                             
      DOUBLE PRECISION Ckn, Cxout, Cphi, Cqtheta                                
      COMMON /imfn1/ ftype                                                      
      COMMON /imfn2/ Ckn, Cxout, Cphi, Cqtheta                                  

      IF (ftype.EQ.1) THEN                                                      
c       this part is for convolution                                            
        x = sqrt(abs(Cxout*Cxout+y*y-2.*Cxout*y*dcos(Cphi)))                    
        imagfn = PSFN(x) * y**Ckn                                               
      ELSE                                                                      
c       this part is for visibility                                             
c       argument is Pi*q*y (not 2*Pi*q*y) to account for the fact that          
c       theta1 is diameter rather than radius (so V is function of              
c       q*theta1, like in IE, '96, MNRAS 279, 1019)                             
        imagfn = Bessel(2.*ASIN(1.0)*Cqtheta*y) * y**Ckn                        
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE MIE(npL,nL,lambda,Ere,Eim,npA,nA,a,nG1,Qabs,Qsca)              
c     =======================================================================       
c     This subroutine calculates Qabs and Qsca for a given diffractive              
c     index Ere, Eim, wavelength lambda and size a. Here, lambda is an              
c     array (1..nL), Ere and Eim are given on this array, a is an array             
c     of sizes (1..nA). Qabs and Qsca are arrays (nG1..nG1+nA,nL), i.e. for         
c     each wavelength lambda, Qabs and Qsca are evaluated for nA different          
c     sizes. The numbering, however, does not start from 1, but rather from         
c     nG1.                                                [Z.I., Aug. 1996]         
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npL, nL, npA, nA, nG1, iL, iA                                     
      DOUBLE PRECISION lambda(npL), Ere(npL), Eim(npL), a(npA),                 
     &       Qabs(npA,npL), Qsca(npA,npL)                                       
      REAL xx, Qex, Qsc, Qback                                                  
      COMPLEX refrel, s1(200), s2(200)                                          

c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c       complex index of refraction                                             
        refrel = cmplx(Ere(iL),Eim(iL))                                         
c       loop over sizes                                                         
        DO iA = 1, nA                                                           
c         size parameter                                                        
          xx=2.0*3.14159265*a(iA)/lambda(iL)                                    
c         if size parameter xx>100 use xx=100 (geometrical optics)              
          IF (xx.GT.100.0) xx = 100.0                                           
c         calculate efficiencies                                                
          CALL bhmie(xx,refrel,2,s1,s2,Qex,Qsc,Qback)                           
c         store the result                                                      
          Qabs(nG1+iA-1,iL) = Qex - Qsc                                         
          Qsca(nG1+iA-1,iL) = Qsc                                               
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      DOUBLE PRECISION FUNCTION PSFN(x)                                         
c     =======================================================================       
c     This function evaluates the point spread function. For psftype.EQ.1           
c     the function is evaluated as a sum of two Gaussians, for psftype.EQ.2         
c     it is provided by user in a file. psftype and all other relevant              
c     parameters come from COMMON /psf/ and are initialized in subroutine           
c     INPUT.                                               [Z.I., Jan. 1997]        
c     =======================================================================       
      IMPLICIT none                                                             
      DOUBLE PRECISION x                                                        
      INTEGER idummy                                                            
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         

      IF (psftype.LT.3) THEN                                                    
        psfn = dexp(-(1.665*x/FWHM1(iLambda))**2.)                              
        IF (psftype.EQ.2)                                                       
     &    psfn = (psfn + kPSF(iLambda) *                                        
     &           dexp(-(1.665*x/FWHM2(iLambda))**2.))/(1.+kPSF(iLambda))        
        ELSE                                                                    
        CALL LinInter(1000,Npsf,xpsf,ypsf,x,idummy,psfn)                        
      ENDIF                                                                     

      RETURN                                                                    
      END                                                                       


      SUBROUTINE setupETA                                                       
c     =======================================================================       
c     This subroutine finds spline coefficients ETAcoef (defined in file            
c     'density.inc') such that normalized density function ETA(Y(iY)) is:           
c     ETAcoef(iY,1)+ETAcoef(iY,2)/Y(iY)+...+ETAcoef(iY,2)/Y(iY)^3                   
c     If spline approximation differs more than maxerr (see below) at the           
c     midpoint, then a straight line is used instead. (In case of wavelength        
c     depend. ETA, use ETAfun where any new dens. laws should be described).        
c     Coefficients ETAcoef are later used in getETAzp to calculate ETAzp.           
c                                                    [ZI, Feb'96; MN,Aug'97]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iY, iC                                                            
      DOUBLE PRECISION coef(npY,4), ETA, maxerr, Ymid, Yinverse(npY),           
     &       ETAaux(npY), ETAmid(npY)                                           

c       generate input function for SPLINE2                                     
        DO iY = 1, nY                                                           
          Yinverse(iY) = 1. / Y(iY)                                             
          ETAaux(iY) = ETA(Y(iY))                                               
          IF (iY.LT.nY) THEN                                                    
            Ymid = dsqrt(Y(iY)*Y(iY+1))                                         
            ETAmid(iY) = ETA(Ymid)                                              
          END IF                                                                
        END DO                                                                  
c       calculate spline coefficients                                           
        CALL SPLINE2(Yinverse,ETAaux,nY,coef)                                   
c       check and fix spline coefficients                                       
        maxerr = 0.1                                                            
c       RDW is initialized in Input   	                                         
        CALL CHKSPLIN(Yinverse,ETAaux,ETAmid,nY,coef,maxerr,RDW)                
c       copy coefficients to the output array ETAcoef                           
        DO iY = 1, nY                                                           
          DO iC = 1, 4                                                          
            ETAcoef(iY,iC) = coef(iY,iC)                                        
          END DO                                                                
        END DO                                                                  

      RETURN                                                                    
      END                                                                       


      DOUBLE PRECISION FUNCTION SizeDist(q,aa,sdtype,a0)                        
c     =======================================================================       
c     This subroutine calculates size distribution n(a) for a=aa. The size          
c     distribution is MRN type n(a)~1/a**q for sdtype.LE.2 and KMH type             
c     n(a)~dexp(-a/a0)/a**q otherwise                                               
c                                                          [Z.I., Aug. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER sdtype                                                            
      DOUBLE PRECISION aa, a0, q                                                

      IF (sdtype.LE.2) THEN                                                     
        SizeDist = 1.0/aa**q                                                    
        ELSE                                                                    
        SizeDist = dexp(-aa/a0)/aa**q                                           
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE Spectral(model,denstyp,nL,Lambda)                              
c     =======================================================================       
c         This subroutine finds the spectral features for spp and zpp files.        
c         It employs SpFeatur().		       		     [MN, Jan'99]                        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER i, iL, model, nL, denstyp, Nchar                                  
      PARAMETER (Nchar=11)                                                      
      DOUBLE PRECISION Charac(Nchar),Spectr(npL),Lambda(npL)                    

c     find features for *.spp file                                              
      IF (denstyp.EQ.0) THEN                                                    
         DO iL = 1, nL                                                          
            Spectr(iL) = ftot(iL,nYok) + ksi*fsR(iL,nYok)                       
       END DO                                                                   
      ELSE                                                                      
         DO iL = 1, nL                                                          
            Spectr(iL) = Spectrum(iL)                                           
       END DO                                                                   
      END IF                                                                    
      CALL SpFeatur(model,nL,Lambda,Spectr,Charac)                              
      DO i = 1, Nchar                                                           
       SpecChar(i,model) = Charac(i)                                            
      END DO                                                                    
c     find features for *.zpp file                                              
      IF (denstyp.EQ.0) THEN                                                    
         DO iL = 1, nL                                                          
            Spectr(iL) = dabs(ftot(iL,1) - fsL(iL,1))                           
         END DO                                                                 
       CALL SpFeatur(model,nL,Lambda,Spectr,Charac)                             
         DO i = 1, Nchar                                                        
            SpecChar(i+11,model) = Charac(i)                                    
         END DO                                                                 
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE SpFeatur(model,nL,Lambda,Spectr,Charac)                        
c     =======================================================================       
c     This subroutine calculates IRAS colors and other spectral quantities          
c     Filters data from Neugebauer et al, 1984, ApJ, 278, L1.                       
c     Procedure described in Bedijn, 1987, A&A, 186, 136.  [Z.I., Mar. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER i, j, iaux, model, nL                                             
      DOUBLE PRECISION  f1(7), f2(7), f3(7), f4(7), w1(7), w2(7), w3(7),        
     &        w4(7), phi(4,7), al(4,7), cl(4), prz(4), tinf(4), flxy(9),        
     &        wav(9), lambda(npL), Spectr(npL), Charac(11), wmid, flx1,         
     &        flx2, faux, B98, B11, rat9818, beta813, beta1422, f12,            
     &        f25, f60, f100, an98, f98c, f11c, an11                            

c     Data for 4 IRAS filters                                                   
c     wavelengths                                                               
       DATA w1/7.55, 8.0, 10.3, 11.5, 13.4, 14.7, 15.5/                         
       DATA w2/16.6, 22.5, 25.6, 26.8, 27.5, 29.3, 31.1/                        
       DATA w3/30.5, 40.1, 40.2, 65.2, 74.2, 83.8, 83.9/                        
       DATA w4/72.7, 95.4, 111.2, 116.6, 137.4, 137.5, 137.6/                   
c     transmittivities                                                          
       DATA f1/0.000, 0.618, 0.940, 0.750, 1.022, 0.906, 0.000/                 
       DATA f2/0.235, 0.939, 0.939, 0.745, 0.847, 0.847, 0.000/                 
       DATA f3/0.000, 0.102, 0.260, 1.026, 0.842, 0.001, 0.000/                 
       DATA f4/0.000, 0.910, 1.000, 0.330, 0.002, 0.001, 0.000/                 
c     ------------------------------------------------------------              
c     initialization                                                            
       DO i = 1, 4                                                              
         tinf(i) = 0.0                                                          
         cl(i) = 0.0                                                            
       END DO                                                                   
       DO j = 1, 7                                                              
         al(1,j) = w1(j)                                                        
         phi(1,j) = f1(j)                                                       
         al(2,j) = w2(j)                                                        
         phi(2,j) = f2(j)                                                       
         al(3,j) = w3(j)                                                        
         phi(3,j) = f3(j)                                                       
         al(4,j) = w4(j)                                                        
         phi(4,j) = f4(j)                                                       
       END DO                                                                   
c     ------------------------------------------------------------------------      
c     first find IRAS colors                                                    
      DO j = 2, nL                                                              
c       middle wavelength                                                       
        wmid = 0.5*(Lambda(j-1)+Lambda(j))                                      
c       interpolate filters for wmid                                            
        CALL PHILAM(wmid,prz,al,phi)                                            
c       convert Spectrum to Flambda                                             
        flx1 = Spectr(j-1) / Lambda(j-1)                                        
        flx2 = Spectr(j) / Lambda(j)                                            
c       add contribution to the integral (index is over filters)                
        DO i = 1, 4                                                             
          tinf(i) = tinf(i) + prz(i)*0.5*(flx2+flx1)*                           
     &                                    (Lambda(j)-Lambda(j-1))               
          cl(i) = cl(i) + prz(i) * (Lambda(j)-Lambda(j-1))                      
        END DO                                                                  
      END DO                                                                    
      DO i = 1, 4                                                               
        tinf(i) = tinf(i) / cl(i)                                               
      END DO                                                                    
c     Spectrum corrected for IRAS filters                                       
      f12 = tinf(1)*12.0                                                        
      f25 = tinf(2)*25.0                                                        
      f60 = tinf(3)*60.0                                                        
      f100 = tinf(4)*100.0                                                      
c     now find  B98, B11, F98/F18, beta 8-13, beta 14-22                        
c     find fluxes at all needed wavelengths (energy increases with index)       
      DATA wav/2.2, 8.0, 9.8, 11.3, 13.0, 14.0, 18.0, 22.0, 0.55/               
      DO j = 1, 9                                                               
        CALL LinInter(npL,nL,lambda,Spectr,wav(j),iaux,faux)                    
        flxy(j) = faux                                                          
      END DO                                                                    
c     the feature strength at 9.8 and 11.4 microns                              
      an98 = log(flxy(5)/flxy(2))/log(wav(5)/wav(2))                            
      f98c = flxy(2)*(wav(3)/wav(2))**an98                                      
      B98 = log(flxy(3)/f98c)                                                   
      an11 = log(flxy(5)/flxy(3))/log(wav(5)/wav(3))                            
      f11c = flxy(3)*(wav(4)/wav(3))**an11                                      
      B11 = log(flxy(4)/f11c)                                                   
c     ratio F9.8/F18                                                            
      rat9818 = flxy(3)/flxy(7)*wav(7)/wav(3)                                   
c     beta 8-13 and beta 14-22 (see Neugebauer)                                 
      beta813 = log(flxy(5)/flxy(2))/log(13.0/8.0) - 1.0                        
      beta1422 = log(flxy(8)/flxy(6))/log(22.0/14.0) - 1.0                      
c     store SpecChar to output array SpecChar                                   
      Charac(1) = B98                                                           
      Charac(2) = B11                                                           
      Charac(3) = rat9818                                                       
      Charac(4) = beta813                                                       
      Charac(5) = beta1422                                                      
      Charac(6) = flxy(9)                                                       
      Charac(7) = flxy(1)                                                       
      Charac(8) = f12                                                           
      Charac(9) = log10(25.0*f25/f12/12.0)                                      
      Charac(10) = log10(60.0*f60/f12/12.0)                                     
      Charac(11) = log10(f100*100.0/f60/60.0)                                   

      RETURN                                                                    
      END                                                                       
                                                                                

      SUBROUTINE PHILAM(Alam,F,Al,Phi)                                          
c     interpolates IRAS filters [f(4)] for a given wavelength alam              
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      DIMENSION Phi(4,7), Al(4,7), F(4)                                         
c     -----------------------------------------------                               
      DO i = 1, 4                                                               
        F(i) = 0.0                                                              
        im = 0                                                                  
        istop = 0                                                               
        DO WHILE (istop.NE.1)                                                   
              im = im + 1                                                       
            IF ( (Alam-Al(i,im))*(Alam-Al(i,im+1)).LE.0) THEN                   
               a = (Phi(i,im+1)-Phi(i,im))                                      
     &           /log10(Al(i,im+1)/Al(i,im))                                    
               b = Phi(i,im) - a*log10(Al(i,im))                                
               F(i) = a*log10(Alam) + b                                         
            IF ( F(i).GT.1.0 ) F(i) = 1.0                                       
              END IF                                                            
            IF (im.EQ.6) istop = 1                                              
        END DO                                                                  
      END DO                                                                    
c     ------------------------------------------------------                        
      RETURN                                                                    
      END                                                                       


      SUBROUTINE Visibili(IntOut)                                               
c     =======================================================================       
c     This subroutine finds visibility functions corresponding to IntOut.           
c     The work horse is subroutine Visi2D, and this subroutine is used to           
c     prepare everything.                                  [Z.I., Jan. 1997]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER i, j, N1, N2                                                      
      DOUBLE PRECISION  IntOut(20,npP+2), Visi(1000), Int1D(npP+2)              

c     generate spatial frequency (q) grid                                       
c     first N1 points up to qtheta1=1.22 (Rayleigh limit for a disk)            
      N1 = 80                                                                   
c     first 2 points manually:                                                  
c     there must be 0!                                                          
      qtheta1(1) = 0.0                                                          
c     make sure the whole envelope is resolved                                  
      qtheta1(2) = 0.5 / bOut(nP+2)                                             
c     and the rest on logarithmic grid up to 1.22                               
      DO i = 1, N1-2                                                            
        qtheta1(i+2) = qtheta1(2) * (1.22/qtheta1(2))**(i*1.0/(N1-2))           
      END DO                                                                    
c     envelope is well sampled, now to be sure that the star will be OK         
c     for small taus add N2 points on a logarithmic grid up to 1.22/p*          
      N2 = 20                                                                   
      DO i = 1, N2                                                              
        qtheta1(N1+i) = 1.22 / bOut(2)**(i*1.0/N2)                              
      END DO                                                                    
      Nvisi = N1 + N2                                                           
c     find visibility wavelength by wavelength                                  
      DO j = 1, NlambdaOut                                                      
        DO i = 1, nP+2                                                          
          Int1D(i) = IntOut(j,i)                                                
          CALL CHKRANGE(dynrange,Int1D(i))                                      
          IF (Int1D(i).LT.dynrange) Int1D(i)=0.0                                
        END DO                                                                  
        CALL Visi2D(npP+2,nP+2,bOut,Int1D,1000,N1+N2,qtheta1,Visi)              
c       copy 1D convolved visibility to Visib                                   
        DO i = 1, N1+N2                                                         
c         check dynamic range                                                   
          CALL CHKRANGE(dynrange,Visi(i))                                       
          Visib(j,i) = Visi(i)                                                  
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE Visi2D(NinMax,Nin,Xin,Yin,Noutmax,Nout,Xout,Yout)              
c     =======================================================================       
c     This subroutine finds the visibility function (the spatial Fourier            
c     transform of the intensity distribution) corresponding to the                 
c     intensity Yin(Xin[i]), i=1,Nin. Visibility, Yout, is evaluated at q           
c     positions (spatial frequency) given in Xout[i], i=1,Nout. Maximum size        
c     of Xin is NinMax, maximum size of Xout is NoutMax. The Bessel function        
c     of the zeroth order is provided separately. The integration is done by        
c     calling subroutine ROMBY (Bessel function is called from IMGFN).              
c     Note:                                                                         
c     The visibility function V(q) for a circularly symmetric intensity             
c     I(x) is:                                                                      
c              V(q) = F(q)/F(0)                                                     
c     where Jo is the Bessel function of the zeroth order, and                      
c              F(q) = Int[Jo(2Pi*q*x)*I(x)*2Pi*x*dx]                                
c     Note that F(0) is nothing more than flux. For more details see                
c     Ivezic & Elitzur, 1996, MNRAS, 279, 1019 and ref. therein.                    
c                                                          [Z.I., Jan. 1997]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER NinMax, Nin, NoutMax, Nout, iq, iXin                              
      DOUBLE PRECISION Xin(NinMax), Yin(NinMax), Xout(NoutMax),                 
     &       Yout(NoutMax),  F(1000), F0, int1, int2, A, B, imagfn              
      INTEGER ftype                                                             
      DOUBLE PRECISION Ckn, Cxout, Cphi, Cqtheta                                
      COMMON /imfn1/ ftype                                                      
      COMMON /imfn2/ Ckn, Cxout, Cphi, Cqtheta                                  
      EXTERNAL imagfn                                                           

c     loop over spatial frequency q (= Xout)                                    
      DO iq = 1, Nout                                                           
        Cqtheta = Xout(iq)                                                      
        F(iq) = 0.0                                                             
c       loop over radial positions                                              
        DO iXin = 1, Nin                                                        
c         find F(q)                                                             
          ftype = 2                                                             
          Ckn = 1.0                                                             
          CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int1)                         
          Ckn = 2.0                                                             
          CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int2)                         
c         contribution from this annulus (lin. approx. for intensity)           
          A = Xin(iXin+1)*Yin(iXin)-Xin(iXin)*Yin(iXin+1)                       
          A = A /(Xin(iXin+1)-Xin(iXin))                                        
          B = (Yin(iXin+1)-Yin(iXin))/(Xin(iXin+1)-Xin(iXin))                   
          F(iq) = F(iq) + A*int1 + B*int2                                       
        END DO                                                                  
      END DO                                                                    
c     flux                                                                      
      F0 = F(1)                                                                 
      DO iq = 1, Nout                                                           
       IF(F0.EQ.0.) THEN                                                        
         Yout(iq) = 0.                                                          
       ELSE                                                                     
         Yout(iq) = dabs(F(iq) / F0)                                            
       END IF                                                                   
      END DO                                                                    

      RETURN                                                                    
      END                                                                       
