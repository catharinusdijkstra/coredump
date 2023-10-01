c =========================================================================     
c     This is are the core subroutines for radiative transfer calculation.      
c                                                             [MN, Mar'99]      
c =========================================================================     
C     Table of Contents                                                         
C                                                                               
C     EMISSION                                                                  
C     GFUN                                                                      
C     FINDTEMP                                                                  
C     GRAYBODY                                                                  
C     INITTEMP                                                                  
C     INSERT                                                                    
C     INTETA                                                                    
C     INVERT                                                                    
C     KINT4                                                                     
C     MATRIX                                                                    
C     NORDLUND                                                                  
C     OCCLTMSG                                                                  
C     PGRID                                                                     
C     PLANCK                                                                    
C     RADTRANSF                                                                 
C     SETGRIDS                                                                  
C     SOLVE                                                                     
C     STAR                                                                      
C     TRAPZD2                                                                   
C     TWOFUN                                                                    
C     WEIGHTS                                                                   
C     YGRID                                                                     
c =========================================================================     
                                                                                

      SUBROUTINE Emission(flag,nG,Td,alpha,abund,Uin,Emiss)                     
c     =======================================================================       
c     This subroutine calculates emission term from the temperature, abund          
c     and alpha arrays for flag=0, and adds U to it for flag=1.                     
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
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER iG, iY, iL, nG, flag                                              
      DOUBLE PRECISION Uin(npL,npY), Td(npG,npY), alpha(npG,npY),               
     &       Emiss(npL,npY), abund(npG,npY), EmiG, xP, Planck                   

c     first initialize Emiss                                                    
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c       loop over radial coordinate                                             
        DO iY = 1, nY                                                           
          Emiss(iL,iY) = 0.0                                                    
        END DO                                                                  
      END DO                                                                    
c     calculate emission term for each component and add it to Emiss            
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c       loop over radial coordinate                                             
        DO iY = 1, nY                                                           
c         loop over grains                                                      
          DO iG = 1, nG                                                         
            xP = 14400.0 / lambda(iL) / Td(iG,iY)                               
            EmiG = abund(iG,iY) * alpha(iG,iY) * Planck(xP)                     
c           add contribution for current grains                                 
            Emiss(iL,iY) = Emiss(iL,iY) + EmiG                                  
          END DO                                                                
c         if needed add Uin                                                     
          IF (flag.EQ.1) THEN                                                   
            Emiss(iL,iY) = Emiss(iL,iY) + Uin(iL,iY)                            
          END IF                                                                
          IF (Emiss(iL,iY).LT.dynrange*dynrange) Emiss(iL,iY) = 0.0             
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       

                                                                                

      DOUBLE PRECISION FUNCTION Gfun(Tt)                                        
c =======================================================================       
c This is an auxiliary function used in determining the dust temperature        
c                                                      [Z.I., Nov. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER iW, nW                                                            
      DOUBLE PRECISION Tt, wav(npL), Eff(npL), ff(npL), xP, gg, f1, f2,         
     &       Planck                                                             
      COMMON /gfunction/ wav, Eff, f1, f2, nW                                   

      DO iW = 1, nW                                                             
        xP = 14400.0 / wav(iW) / Tt                                             
        ff(iW) = Eff(iW) * Planck(xP) / wav(iW)                                 
      END DO                                                                    
      CALL Simpson(npL,1,nW,wav,ff,gg)                                          
      Gfun = f1 - f2 * gg * Tt**4.0                                             

      RETURN                                                                    
      END                                                                       

                                                                                
      SUBROUTINE FindTemp(geom,Utot,nG,Td,alpha)                                
c     =======================================================================       
c     This subroutine finds new temperature and alpha arrays from Utot.             
c     Temperature is obtained by solving:                                           
c       f1(y) - f2(y)*g(T) = 0                                                      
c     where                                                                         
c       f1(y) = Int(Qabs*Utot*dlambda)                                              
c       f2(y) = alpha(y=1)*y**2/Tsub**4                                             
c        g(T) = Td**4*Int(Qabs*Planck(Td)*dlambda)                                  
c     alpha(y)= Int(Qabs*Utot*dlambda)/Int(Qabs*Planck(Td)*dlambda)                 
c     The equation is solved by Ridder's method (subroutine from Num.Rec)           
c                                                          [Z.I., Oct. 1996]        
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
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER iG, iY, iL,nG, nW, succes, geom                                   
      DOUBLE PRECISION Utot(npL,npY), Td(npG,npY), alpha(npG,npY),              
     &       fnum(npL), fdenum(npL), num, denum, xP, Tin, aa, Planck,           
     &       wav(npL), f1, f2, gfun, x1, x2,  Eff(npL), Zriddr                  
      EXTERNAL gfun                                                             
      COMMON /gfunction/ wav, Eff, f1, f2, nW                                   

      nW = nL                                                                   
c     loop over grains                                                          
      DO iG = 1, nG                                                             
c       evaluate alpha(1)                                                       
        DO iL = 1, nL                                                           
          fnum(iL) = SigmaA(iG,iL) * Utot(iL,1) / lambda(iL)                    
          xP = 14400.0 / lambda(iL) / Tsub(iG)                                  
          fdenum(iL) = SigmaA(iG,iL) * Planck(xP) / lambda(iL)                  
          wav(iL) = lambda(iL)                                                  
          Eff(iL) = SigmaA(iG,iL)                                               
        END DO                                                                  
        CALL Simpson(npL,1,nL,lambda,fnum,num)                                  
        CALL Simpson(npL,1,nL,lambda,fdenum,denum)                              
        alpha(iG,1) = num /  denum                                              
        aa = alpha(iG,1) / Tsub(iG)**4.0                                        
        Td(iG,1) = Tsub(iG)                                                     
c       loop over radial positions (solving f1-f2*g(T)=0)                       
        DO iY = 2, nY                                                           
c         calculate f1 and f2                                                   
          DO iL = 1, nL                                                         
            fnum(iL) = SigmaA(iG,iL) * Utot(iL,iY) / lambda(iL)                 
          END DO                                                                
          CALL Simpson(npL,1,nL,lambda,fnum,f1)                                 
          IF (geom.EQ.0) THEN                                                   
c           This is for slab geometry ...                                       
            f2 = aa                                                             
          ELSE                                                                  
c          ... and for spherical                                                
             f2 = aa * Y(iY)**2.0                                               
          END IF                                                                
c         estimate range for temperature                                        
          x1 = Td(iG,iY)                                                        
          x2 = 1.01*Td(iG,iY)                                                   
c         make sure that correct solution is bracketted                         
          CALL Zbrac(gfun,x1,x2,100,succes)                                     
          IF (succes.NE.1) THEN                                                 
            CALL Msg(4)                                                         
            ELSE                                                                
c           save the old value of Td                                            
            Tin = Td(iG,iY)                                                     
c           get new temperature                                                 
            Td(iG,iY) = Zriddr(gfun,x1,x2,500,accConv)                          
c           calculate alpha                                                     
            DO iL = 1, nL                                                       
              xP = 14400.0 / lambda(iL) / Td(iG,iY)                             
              fdenum(iL) = SigmaA(iG,iL) * Planck(xP) / lambda(iL)              
            END DO                                                              
            CALL Simpson(npL,1,nL,lambda,fdenum,denum)                          
            alpha(iG,iY) = f1 /  denum                                          
          END IF                                                                
        END DO                                                                  
      END DO                                                                    

999   RETURN                                                                    
      END                                                                       


      SUBROUTINE GrayBody(albedo,TAUgbTot,Ugb,fgb)                              
c     =======================================================================       
c     This subroutine solves the gray body problem for albedo=1 (or                 
c     equivalently pure scattering) and scattering with absorption (but no          
c     emission) for albedo<1, in a spherically symmetric envelope. Total            
c     optical depth is TAUtot, and density law is specified elsewhere.              
c     This subroutine was designed to be a part of Dusty and to use already         
c     existing subroutines as much as possible, so some parts might seem to         
c     be a little awkward.                           [ZI,Jul'96;MN,Sep'97]          
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
      INTEGER iPGB, iY, nLstore, error                                          
      DOUBLE PRECISION ETAzp(npP,npY), mat0(npL,npY,npY),                       
     &       mat1(npL,npY,npY), Em(npL,npY), Ugb(npY), fgb(npY),                
     &       Dummy1(npL,npP,npY), Dummy2(npL,npP,npY), Dummy3(npL,npY),         
     &       albedo, pGB, TAUgbTot, TAUstore                                    
c ----------------------------------------------------------------------        
c     Values needed in this subroutine only                                     
      pGB = 0.0                                                                 
      iPGB = 0                                                                  
      nLstore = nL                                                              
      nL = 1                                                                    
      TAUstore = TAUtot(1)                                                      
      TAUtot(1) = TAUgbTot                                                      
c     generate spline coefficients for ETA                                      
      CALL setupETA                                                             
c     evaluate ETAzp                                                            
      CALL getETAzp(ETAzp)                                                      
c     generate some temporary arrays                                            
      DO iY = 1, nY                                                             
        Us(1,iY) = dexp(-ETAzp(1,iY)*TAUgbTot)                                  
        fs(1,iY) = Us(1,iY)                                                     
        Em(1,iY) = 0.0                                                          
        fde(1,iY) = 0.0                                                         
        omega(1,iY) = albedo                                                    
      END DO                                                                    
c     find radiative transfer matrices                                          
      CALL Matrix(ETAzp,pGB,iPGB,mat0,mat1,Dummy1,Dummy2)                       
c     solve for Utot                                                            
      CALL Invert(nY,nL,mat0,Us,Em,omega,Utot,Dummy3,error)                     
c     calculate flux, ftot                                                      
      CALL Multiply(1,npY,nY,npL,nL,mat1,Utot,omega,1,fs,ftot)                  
c     store to the output arrays                                                
      DO iY = 1, nY                                                             
        Ugb(iY) = Utot(1,iY)                                                    
        fgb(iY) = ftot(1,iY)                                                    
      END DO                                                                    
      nL = nLstore                                                              
      TAUtot(1) = TAUstore                                                      

      RETURN                                                                    
      END                                                                       


      SUBROUTINE InitTemp(ETAzp,nG,alpha)                                       
c     =======================================================================       
c     This subroutine initializes the first approximations for temperature          
c     and alpha arrays. It uses approximations given by eqs. B5 and B7 from         
c     IE96. Alpha(y) is defined as the ratio of Qabs averaged with the total        
c     energy density at y and Qabs Planck averaged with T(y). That is, Psi          
c     is alpha(y=1) and (T/Tsub)^4 = alpha(y) / alpha(1).  [Z.I., Jul. 1996]        
c                                                                                   
c     *** This version works for single grains ONLY ***                             
c                                                                                   
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
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      INTEGER iL, iY, iG, nG, istop, iW                                         
      DOUBLE PRECISION ETAzp(npP,npY), xP, Planck, resaux, aux(npY),            
     &       faux(npL), fstar(npY), Qfstar(npY),  fstarQ(npY),                  
     &       alpha(npG,npY), oldalpha(npY), Sigext(npL), QP(npY),               
     &       QF(npY),EtaQF, ETA, delta                                          

c     Lambda integral of fs=f_e*exp(-TAU) -> fstar                              
      DO iY = 1, nY                                                             
c       generate auxiliary function for lambda integration:                     
        DO iL = 1, nL                                                           
          faux(iL) = fs(iL,iY) / lambda (iL)                                    
        END DO                                                                  
        CALL Simpson(npL,1,nL,lambda,faux,resaux)                               
        fstar(iY) = resaux                                                      
      END DO                                                                    
c     loop over grains                                                          
      DO iG = 1, nG                                                             
c       first approximation for temperature                                     
        DO iY = 1, nY                                                           
          Td(iG,iY) = Tsub(iG) / Y(iY)**0.4                                     
          alpha(iG,iY) = 1.0                                                    
        END DO                                                                  
c       generate the extinction cross-section Sigext                            
        DO iL = 1, nL                                                           
          Sigext(iL) = SigmaA(iG,iL)+SigmaS(iG,iL)                              
        END DO                                                                  
c       Lambda integral of Sigext*fs -> Qfstar in IE97                          
        DO iY = 1, nY                                                           
c         generate auxiliary function for lambda integration:                   
          DO iL = 1, nL                                                         
            faux(iL) = Sigext(iL) * fs(iL,iY) / lambda (iL)                     
          END DO                                                                
          CALL Simpson(npL,1,nL,lambda,faux,resaux)                             
          Qfstar(iY) = resaux                                                   
        END DO                                                                  
c       iterate until Td and Psi converge (i.e. until alpha does)               
        istop = 0                                                               
        DO WHILE (istop.NE.1)                                                   
          istop = 1                                                             
c         find Planck average of Sigext and calculate QF (eq. B5)               
          DO iY = 1, nY                                                         
c           generate auxiliary function for lambda integration:                 
            DO iL = 1, nL                                                       
              xP = 14400.0 / Td(iG,iY) / lambda(iL)                             
              faux(iL) = Sigext(iL) * Planck(xP) / lambda (iL)                  
            END DO                                                              
            CALL Simpson(npL,1,nL,lambda,faux,resaux)                           
c           Planck average of Sigext                                            
            QP(iY) = resaux                                                     
c           calculate intetgral next to 1/y2 in eq. B7 -> fstarQ                
            DO iL = 1, nL                                                       
              faux(iL) = fs(iL,iY)*(1.-SigmaA(iG,iL)/QP(iY))/lambda(iL)         
            END DO                                                              
            CALL Simpson(npL,1,nL,lambda,faux,resaux)                           
            fstarQ(iY) = resaux                                                 
c           calculate QF (eq. B5)                                               
            QF(iY) = Qfstar(iY) + QP(iY) * (1.0-fstar(iY))                      
c           store current alpha                                                 
            oldalpha(iY) = alpha(iG,iY)                                         
          END DO                                                                
c         for each y calculate new alpha from eq. B7                            
c         two y loops are needed because of integral QF*ETA                     
          DO iY = 1, nY                                                         
c           first evaluate integral of QF with ETA...                           
            DO iW = iY, nY                                                      
              aux(iW) = QF(iW) * ETA(Y(iW)) / Y(iW)**2.0                        
            END DO                                                              
            CALL Simpson(npY,iY,nY,Y,aux,EtaQF)                                 
c           ... and alpha                                                       
            alpha(iG,iY) = 3.*TAUfid*EtaQF + (1-fstarQ(iY))                     
c           calculate temperature                                               
            Td(iG,iY) = Tsub(iG) * (alpha(iG,iY)/alpha(iG,1))**0.25             
            Td(iG,iY) = Td(iG,iY) / dsqrt(Y(iY))                                
            delta = DABS((alpha(iG,iY)-oldalpha(iY))/alpha(iG,iY))              
            IF (delta.GT.accConv) istop = 0                                     
          END DO                                                                
        END DO                                                                  
c       end of iterations                                                       
      END DO                                                                    
c     end of loop over grains (iG)                                              

      RETURN                                                                    
      END                                                                       

                                                                                
      SUBROUTINE INSERT(pstar,iP,iPstar)                                        
c     =======================================================================       
c     This subroutine inserts two rays if the stellar disk is finite. The           
c     first ray, corresponding to the disk edge, has index iPstar, and the          
c     following ray with 1.0001 greater impact parameter has index iPstar+1.        
c     The only exception is if pstar>0.9999 when only a single ray is               
c     inserted.                                            [Z.I., Feb. 1996]        
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
      INTEGER iP, iPstar                                                        
      DOUBLE PRECISION pstar                                                    

      IF (pstar.GE.0.999999) THEN                                               
        pstar = 0.999999                                                        
        iP = nPcav + 1                                                          
        iPstar = nPcav + 1                                                      
        P(iP) = pstar                                                           
        ELSE                                                                    
        IF (pstar.GE.P(nPcav)) THEN                                             
          IF (pstar.EQ.P(nPcav)) pstar = 1.00001*P(nPcav)                       
          P(nPcav+1) = pstar                                                    
          iPstar = nPcav + 1                                                    
          P(nPcav+2) = 1.00001*pstar                                            
          iP = nPcav + 2                                                        
        ELSE                                                                    
          iPstar = 0                                                            
          iP = 0                                                                
          DO WHILE (iPstar.EQ.0)                                                
            iP = iP + 1                                                         
            IF (P(iP).GT.pstar.AND.iPstar.EQ.0) iPstar = iP                     
          END DO                                                                
          DO iP = 1, nPcav-iPstar+1                                             
            P(nPcav+3-iP) = P(nPcav+1-iP)                                       
          END DO                                                                
          IF (pstar.EQ.P(iPstar-1)) pstar = 1.00001*P(iPstar-1)                 
          P(iPstar) = pstar                                                     
          P(iPstar+1) = 1.00001*pstar                                           
          iP = nPcav + 1                                                        
        END IF                                                                  
      END IF                                                                    

      RETURN                                                                    
      END                                                                       

                                                                                

      DOUBLE PRECISION FUNCTION IntETA(p,iW1,w1,w)                              
c =======================================================================       
c This function calculates the integral over the normalized dens. prof.         
c along the line of sight with impact parameter p and between the points        
c corresponding to y=w1 and y=w. The method used is spline approximation        
c for normalized density distribution ETA and subsequent integration            
c performed analytically by MAPLE (these results are given through              
c soubroutine Maple3).                         [ZI,Feb'96,MN,Aug'97]            
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
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
      INTEGER iW1, iC                                                           
      DOUBLE PRECISION  p, w1, w, aux(4), z, z1, aux1(4)                        

      z = dsqrt(w*w-p*p)                                                        
      z1 = dsqrt(w1*w1-p*p)                                                     
c     integrals calculated by MAPLE                                             
      CALL Maple3(w,z,p,aux)                                                    
      CALL Maple3(w1,z1,p,aux1)                                                 
      DO iC = 1, 4                                                              
        aux(iC) = aux(iC) - aux1(iC)                                            
      END DO                                                                    
      IntETA = 0.0                                                              
      DO iC = 1, 4                                                              
        IntETA = IntETA + ETAcoef(iW1,iC) * aux(iC)                             
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE INVERT(nY,nL,mat,Us,Em,omega,Utot,Uold,error)                  
c     =======================================================================       
c     This subroutine solves the linear system                                      
c     [Utot] = [Us+Em] + [mat0]*[omega*Utot] by calling LINSYS subroutine.          
c                                                          [Z.I., Nov. 1995]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER nY, iL,nL,iY, iYaux, Kronecker, error                             
      DOUBLE PRECISION mat(npL,npY,npY), Us(npL,npY), Em(npL,npY),              
     &       Utot(npL,npY), Uold(npL,npY),B(npY), A(npY,npY), X(npY),           
     &       omega(npL,npY)                                                     

      error = 0                                                                 
c     first copy Utot to Uold                                                   
      DO iL = 1, nL                                                             
        DO iY = 1, nY                                                           
          Uold(iL,iY) = Utot(iL,iY)                                             
        END DO                                                                  
      END DO                                                                    
c     calculate new energy density                                              
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c       generate the vector of free coefficients, B, and matrix A               
        DO iY = 1, nY                                                           
          B(iY) = Us(iL,iY)                                                     
          DO iYaux = 1, nY                                                      
            Kronecker = 0                                                       
            IF (iY.EQ.iYaux) Kronecker = 1                                      
            B(iY) = B(iY) +                                                     
     &                (1.-omega(iL,iYaux))*Em(iL,iYaux)*mat(iL,iY,iYaux)        
            A(iY,iYaux) = Kronecker - omega(iL,iYaux) * mat(iL,iY,iYaux)        
          END DO                                                                
        END DO                                                                  
                                                                                
c       solve the system                                                        
        CALL LINSYS(nY,A,B,X,error)                                             
        IF(error.NE.0) THEN                                                     
         CALL MSG(20)                                                           
         iERROR = iERROR + 1                                                    
         RETURN                                                                 
        END IF                                                                  
c       store the result                                                        
        DO iY = 1, nY                                                           
          IF (X(iY).GE.dynrange*dynrange) THEN                                  
            Utot(iL,iY) = X(iY)                                                 
          ELSE                                                                  
            Utot(iL,iY) = 0.0                                                   
          END IF                                                                
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE Kint4(TAUaux,iP,iL,nZ,K1p,K2p,K3p,K4p,K1m,K2m,K3m,K4m)         
c     =======================================================================       
c     For given wavelength (iL) and impact parameter (iP), this subroutine          
c     calculates integrals Knp and Knm defined as:                                  
c                   Knp(iZ)=INT[PHIn(tz)*exp(tz)*dtz]                               
c     from tz1=TAUaux(iL,iP,iZ) to tz2=TAUaux(iL,iP,iZ+1) and analogously for       
c     Km with exp(tz) replaced by exp(-tz). Function PHIn is defined as             
c     x**(n-1)/Yloc^2, where Yloc is the local radius corresponding to tz,          
c     and x measures relative radial tau: x = (rt - tL)/(tR-tL). Here rt is         
c     the radial optical depth corresponding to tz and tL and tR are radial         
c     optical depths at the boundaries of the integration interval:                 
c     tL = TAUaux(iL,1,iZ) = rt(iZ) and tR = TAUaux(iL,1,iZ+1) = rt(iZ+1).          
c     Integration is performed in z space by Romberg integration implemented        
c     in subroutine ROMBERG2 (slightly changed version of 'qromb' from Num.         
c     Recipes).                                     [ZI,Feb'96;MN,Sep'97]           
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
      INTEGER iP, iL, nZ, iZ, iW1, iLaux, iC                                    
      DOUBLE PRECISION TAUaux(npL,npP,npY), K1p(npY),K2p(npY),K3p(npY),         
     &    K4p(npY), K1m(npY), K2m(npY), K3m(npY), K4m(npY), result(8),          
     &    Kaux(8), deltrton(4),tRL,paux, w1,w2,wL, delTAUzp, z1, z2             
      COMMON /phi2/ paux, w1, wL, delTAUzp, iLaux, iW1                          

      paux = P(iP)                                                              
      iLaux = iL                                                                
c     iLaux is needed to avoid compiler errors since it is in COMMON            
c     /phi2/ (here and in 'TWOfun'), while iL is transferred as a               
c     argument loop over positions on the line of sight                         
      DO iZ = 1, nZ-1                                                           
c       index for the local radial position (left boundary)                     
        iW1 = iYfirst(iP) + iZ - 1                                              
c       radii at the boundaries                                                 
        wL = Y(iW1)                                                             
        IF (iZ.EQ.1) THEN                                                       
          if (paux.GT.1.0) then                                                 
            w1 = paux                                                           
            else                                                                
            w1 = 1.0                                                            
          end if                                                                
          ELSE                                                                  
          w1 = Y(iW1)                                                           
        END IF                                                                  
        w2 = Y(iW1+1)                                                           
        z1 = dsqrt(DABS(w1*w1 - paux*paux))                                     
        z2 = dsqrt(w2*w2 - paux*paux)                                           
c       radial tau-difference at the bound., scaled to tot. opt. depth          
        tRL = TAUaux(iL,1,iW1+1)-TAUaux(iL,1,iW1)                               
c       auxiliary quantity aux/tRL**(n-1)                                       
        deltrton(1) = TAUtot(iL)                                                
        DO iC = 1, 3                                                            
          deltrton(iC+1) = deltrton(iC)/tRL                                     
        END DO                                                                  
c       delTAUzp is needed in PHIn fun's                                        
        delTAUzp = TAUaux(iL,iP,iZ+1)-TAUaux(iL,iP,iZ)                          
c       integrate this step for all 8 cases                                     
        CALL ROMBERG2(z1,z2,result)                                             
c       generate output values                                                  
        DO iC = 1, 4                                                            
           Kaux(iC) = result(iC) * deltrton(iC)                                 
           Kaux(iC+4) = result(iC+4) * deltrton(iC)                             
        END DO                                                                  
        K1m(iZ) = Kaux(1)                                                       
        K2m(iZ) = Kaux(2)                                                       
        K3m(iZ) = Kaux(3)                                                       
        K4m(iZ) = Kaux(4)                                                       
        K1p(iZ+1) = Kaux(5)                                                     
        K2p(iZ+1) = Kaux(6)                                                     
        K3p(iZ+1) = Kaux(7)                                                     
        K4p(iZ+1) = Kaux(8)                                                     
      END DO                                                                    
c     set undefined elements to 0                                               
      K1m(nZ) = 0.0                                                             
      K2m(nZ) = 0.0                                                             
      K3m(nZ) = 0.0                                                             
      K4m(nZ) = 0.0                                                             
      K1p(1) = 0.0                                                              
      K2p(1) = 0.0                                                              
      K3p(1) = 0.0                                                              
      K4p(1) = 0.0                                                              

      RETURN                                                                    
      END                                                                       


      SUBROUTINE matrix(ETAzp,pstar,iPstar,m0,m1,mifront,miback)                
c     =======================================================================       
c     This subroutine evaluates radiative transfer matrix for spherically           
c     symmetric envelope. Here m is the order of moment (0 for energy dens.,        
c     1 for flux, 2 for pressure etc.), ETAzp is array of optical depths            
c     along the line of sight and mat is radiative transfer matrix.                 
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER m                                                                 
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
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
      INTEGER iP, iL, nZ(npP), iZ, iW,  iY, flag, im, iPstar, jZ, error         
      DOUBLE PRECISION ETAzp(npP,npY),TAUaux(npL,npP,npY), haux(npP),H,         
     &       m0(npL,npY,npY), m1(npL,npY,npY), pstar, addplus,addminus,         
     &       Tplus(npP,npY,npY),Tminus(npP,npY,npY), xN(npP),yN(npP),           
     &       result1,result2,resaux,TAUr(npY), wm(npY),wmT(npY),wp(npY),        
     &       alpha(npY,npY),beta(npY,npY),gamma(npY,npY),delta(npY,npY),        
     &       wgmatp(npY,npY), wgmatm(npY,npY), mifront(npL,npP,npY),            
     &       miback(npL,npP,npY), fact, faux                                    

      error = 0                                                                 
c     generate auxiliary arrays haux & nZ                                       
      DO iP = 1, nP                                                             
c       parameter alowing for a line of sight terminating on the star           
c       H(x1,x2) is the step function.                                          
        haux(iP) = H(P(iP),pstar)                                               
c       upper limit for the counter of z position                               
        nZ(iP) = nY + 1 - iYfirst(iP)                                           
      END DO                                                                    
c     Using the local array TAUaux to avoid multiple calculations of the        
c     product                                                                   
      DO iL = 1, nL                                                             
       DO iP = 1, nP                                                            
        DO iY = 1, nY                                                           
           TAUaux(iL,iP,iY) = ETAzp(iP,iY)*TAUtot(iL)                           
        END DO                                                                  
       END DO                                                                   
      END DO                                                                    
c     -- evaluate matrix elements --                                            
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c       radial optical depths                                                   
        DO iY = 1, nY                                                           
          TAUr(iY) = ETAzp(1,iY)*TAUtot(iL)                                     
        END DO                                                                  
c       auxiliary arrays for given TAUr                                         
        CALL MYSPLINE(TAUr,nY,alpha,beta,gamma,delta)                           
c       loop over impact parameters                                             
        DO iP = 1, nP-1                                                         
c         set T-s to 0                                                          
          DO iY = 1, nY                                                         
            DO iW = 1, nY                                                       
              Tplus(iP,iY,iW) = 0.0                                             
              Tminus(iP,iY,iW) = 0.0                                            
            END DO                                                              
          END DO                                                                
c         generate weights matrices                                             
          CALL WEIGHTS(TAUaux,iP,iL,nZ(iP),nY,alpha,beta,gamma,delta,           
     &                 wgmatp,wgmatm)                                           
c         first position on the line of sight                                   
          IF (YPequal(iP).EQ.1) THEN                                            
            iY = iYfirst(iP)                                                    
            DO iW = 1, nY                                                       
c             cummulative weights for parts II & III                            
              wmT(iW) = 0.0                                                     
              DO jZ = 1, nZ(iP)-1                                               
                fact = dexp(-TAUaux(iL,iP,jZ))                                  
                wmT(iW) = wmT(iW) + fact * wgmatm(jZ,iW)                        
              END DO                                                            
              Tplus(iP,iY,iW) = Tplus(iP,iY,iW) + haux(iP)*wmT(iW)              
              Tminus(iP,iY,iW) = Tminus(iP,iY,iW) + wmT(iW)                     
            END DO                                                              
          END IF                                                                
c         loop over positions on the line of sight                              
          DO iZ = 2, nZ(iP)                                                     
c           increase index for radial position                                  
            iY = iYfirst(iP) + iZ - 1                                           
c           generate weights for this position                                  
            DO iW = 1, nY                                                       
              wp(iW) = 0.0                                                      
              wm(iW) = 0.0                                                      
              wmT(iW) = 0.0                                                     
c             part I                                                            
              DO jZ = 2, iZ                                                     
                fact = dexp(TAUaux(iL,iP,jZ)-TAUaux(iL,iP,iZ))                  
                wp(iW) = wp(iW) + fact * wgmatp(jZ,iW)                          
              END DO                                                            
c             part II & III                                                     
              DO jZ = 1, nZ(iP)-1                                               
                fact = dexp(-(TAUaux(iL,iP,iZ)+TAUaux(iL,iP,jZ)))               
                wmT(iW) = wmT(iW) + fact * wgmatm(jZ,iW)                        
              END DO                                                            
c             part IV                                                           
              IF (iZ.LT.nZ(iP)) THEN                                            
                DO jZ = iZ, nZ(iP)-1                                            
                  fact = dexp(-(TAUaux(iL,iP,jZ)-TAUaux(iL,iP,iZ)))             
                  wm(iW) = wm(iW) + fact * wgmatm(jZ,iW)                        
                END DO                                                          
              ELSE                                                              
                wm(iW) = 0.0                                                    
              END IF                                                            
c             add contribution from this step                                   
              addplus = wp(iW) + haux(iP)*wmT(iW)                               
              Tplus(iP,iY,iW) = Tplus(iP,iY,iW) + addplus                       
              addminus = wm(iW)                                                 
              Tminus(iP,iY,iW) = Tminus(iP,iY,iW) + addminus                    
            END DO                                                              
c         end of loop over iZ                                                   
          END DO                                                                
c       end of the impact parameter loop, iP                                    
        END DO                                                                  
c       add points on the edge                                                  
        DO iW = 1, nY                                                           
          Tplus(nP,nY,iW) = 0.0                                                 
          Tminus(nP,nY,iW) = 0.0                                                
        END DO                                                                  
c     ============================                                              
c       find mat(iL,iY,iW) -> angular (mu) integration                          
c       loop over moments (without calculation of rad.pressure)                 
        DO im = 1, 2                                                            
          m = im - 1                                                            
c         loop over radial positions                                            
          DO iY = 1, nY                                                         
c           generate mu arrray                                                  
            DO iP = 1, Plast(iY)                                                
              xN(iP) = sqrt(1.0-(P(iP)/Y(iY)*P(iP)/Y(iY)))                      
            END DO                                                              
c           loop over local (radial) positions                                  
            DO iW = 1, nY                                                       
c             generate intensity array for NORDLUND                             
              DO iP = 1, Plast(iY)                                              
c               'faux' is a representation of (-1)**m                           
                 faux = 1.0 - 2.0*MOD(m,2)                                      
                 yN(iP) = Tplus(iP,iY,iW) + faux*Tminus(iP,iY,iW)               
c               store matrix elements to evaluate intensity (*1/4Pi)            
                IF (im.EQ.1.AND.iY.EQ.nY) THEN                                  
                  mifront(iL,iP,iW) = 0.0795775 * Tplus(iP,iY,iW)               
                  miback(iL,iP,iW) = 0.0795775 * Tminus(iP,iY,iW)               
                END IF                                                          
              END DO                                                            
c             angular integration inside cavity                                 
              IF (pstar.GT.0.0) THEN                                            
                CALL NORDLUND(0,xN,yN,1,iPstar,m,resaux,error)                  
                IF (error.NE.0) GOTO 999                                        
                IF (nPcav.GT.iPstar)                                            
     &          CALL NORDLUND(0,xN,yN,iPstar+1,nPcav+1,m,result1,error)         
                IF (error.NE.0) GOTO 999                                        
                result1 = result1 + resaux                                      
              ELSE                                                              
                CALL NORDLUND(0,xN,yN,1,nPcav+1,m,result1,error)                
                IF (error.NE.0) GOTO 999                                        
              END IF                                                            
c             flag for analytic integration outside cavity                      
              IF (iY.GT.6) THEN                                                 
                flag = 1                                                        
              ELSE                                                              
                flag = 0                                                        
              ENDIF                                                             
c             angular integration outside cavity                                
              IF (iY.GT.1) THEN                                                 
                CALL NORDLUND(flag,xN,yN,nPcav+1,Plast(iY),                     
     &                                      m,result2,error)                    
                IF (error.NE.0) GOTO 999                                        
              ELSE                                                              
                result2 = 0.0                                                   
              END IF                                                            
c             store current matrix element                                      
              IF (m.EQ.0)                                                       
     &        m0(iL,iY,iW) = 0.5*Y(iY)*Y(iY)*(result1 + result2)                
              IF (m.EQ.1)                                                       
     &        m1(iL,iY,iW) = 0.5*Y(iY)*Y(iY)*(result1 + result2)                
            END DO                                                              
          END DO                                                                
        END DO                                                                  
c     =============================                                             
c     end of loop over wavelengths                                              
      END DO                                                                    
c     save Y and P grids to Yok and Pok, they are needed for analysis           
c     in cases when requirement for finer grids cannot be satisfied and         
c     previous solution is used for output                                      
      nYok = nY                                                                 
      DO iY = 1, nY                                                             
        Yok(iY) = Y(iY)                                                         
      END DO                                                                    
      nPok = nP                                                                 
      DO iP = 1, nP                                                             
        Pok(iP) = P(iP)                                                         
      END DO                                                                    

999   RETURN                                                                    
      END                                                                       


      SUBROUTINE NORDLUND(flag,x,y,N1,N2,m,intydx,error)                        
c     =======================================================================       
c     This subroutine calculates integral I(x**m*y(x)*dx). Both y and x are         
c     1D arrays, y(i), x(i) with i=1,nM (nM comes from 'paramslb.inc'). Lower       
c     and upper integration limits are x(N1) and x(N2), respectively. The           
c     method used is approximation of y(x) by a piecewise cubic spline (see         
c     Nordlund: Spherical Transfer with Single-Ray Approximation, in                
c     'Methods in radiative transfer', ed. W. Kalkofen, Cambridge University        
c     Press, 1984). The resulting integral is sum of y(i)*w(i), i=N1,N2.            
c     Weights w(i) are determined from the array x(i). Here w(i) = wSimp(i)         
c     + wCorr(i). To improve accuracy of angular integration in radiative           
c     transfer, if flag=1 the contribution of the last Nanal (specified             
c     below) points is calculated in subroutine ANALINT which fits a                
c     function of the form: y = P(x) + d/sqrt(1-x*x), where P(x) is the             
c     polynomial of order Nanal-1, to these points and evaluates integral           
c     analytically.                                        [Z.I., Nov. 1995]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER i, flag, N1, N2, N2n, Nanal, m, first, error                      
      DOUBLE PRECISION x(npP), y(npP), wSimp(npP), wCorr(npP), wC1, wC2,        
     &        am, intydx, xaux(4), yaux(4), aux                                 

      error = 0                                                                 
c     parameter 'first' selects choice for derivatives at boundary points.          
c     For first.EQ.0 derivatives are 0 and first*(y2-y1)/(x2-x1) otherwise.         
c     first=1 works much better for functions encountered here.                     
      first = 1                                                                 
c     number of points for analytic integration                                     
      Nanal = 4                                                                 
c     do not include points for analytic integration                                
      IF (flag.EQ.1.AND.N2.GT.N1+Nanal) THEN                                    
        N2n = N2 - Nanal + 1                                                    
      ELSE                                                                      
        N2n = N2                                                                
      END IF                                                                    
c     set integral to 0 and accumulate result in the loop                           
      intydx = 0.0                                                              
c     generate weighting factors, w(i), and integrate in the same loop              
      DO i = N1, N2n                                                            
c       first usual Simpson factors (Nordlund, eq. III-14, first term)...             
        IF (i.NE.N1.AND.i.NE.N2n) THEN                                          
          wSimp(i) = 0.5 * (x(i+1)-x(i-1))                                      
        ELSE                                                                    
          IF (i.eq.N1) wSimp(i) = 0.5 * (x(N1+1)-x(N1))                         
          IF (i.eq.N2n) wSimp(i) = 0.5 * (x(N2n)-x(N2n-1))                      
        END IF                                                                  
c       ... and then correction term for cubic spline (Nordlund, eq. III-14,          
c       second term and eq. III-16) (wC1 and wC2 are auxiliary quantities)            
        IF (i.GT.N1+1) THEN                                                     
          wC1 = x(i) - 2.0*x(i-1) + x(i-2)                                      
        ELSE                                                                    
          IF (i.EQ.N1) wC1 =  first * (x(N1+1) - x(N1))                         
          IF (i.EQ.N1+1) wC1 = first * (x(N1) - x(N1+1))                        
        ENDIF                                                                   
        IF (i.LE.(N2n-2)) THEN                                                  
          wC2 = x(i+2) - 2.0*x(i+1) + x(i)                                      
        ELSE                                                                    
          IF (i.EQ.(N2n-1)) wC2 = first * (x(N2n-1) - x(N2n))                   
          IF (i.EQ.N2n) wC2 = first * (x(N2n) - x(N2n-1))                       
        ENDIF                                                                   
        wCorr(i) = (wC1 - wC2) / 12.                                            
c       add contribution to the integral                                        
        IF (m.EQ.0) THEN                                                        
          intydx = intydx + y(i) * (wSimp(i) + wCorr(i))                        
         ELSE IF(m.EQ.1) THEN                                                   
            intydx = intydx + x(i)*y(i)*(wSimp(i) + wCorr(i))                   
         ELSE IF(m.EQ.2) THEN                                                   
            intydx = intydx + x(i)*x(i)*y(i)*(wSimp(i) + wCorr(i))              
        END IF                                                                  
      END DO                                                                    
c     change the sign (x [i.e. mu] array is in descending order!!!)             
      intydx = -intydx                                                          
c     if flag=1 use analytic approximation for the last Nanal points            
      IF (flag.EQ.1.AND.N2n.GT.N1+Nanal) THEN                                   
c       generate auxiliary arrays for ANALINT                                   
        DO i=1,Nanal                                                            
          xaux(i) = x(N2n+Nanal-i)                                              
          yaux(i) = y(N2n+Nanal-i)                                              
        END DO                                                                  
c       calculate the contribution of the last Nanal points                       
c       produce REAL copy of m                                                  
        am = 1.0*(m)                                                            
        CALL ANALINT(Nanal,xaux,yaux,am,aux,error)                              
        IF(error.NE.0) RETURN                                                   
c       add the contribution of the last Nanal points                           
        intydx = intydx + aux                                                   
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE OccltMSG                                                       
c     =======================================================================       
c         Prints a message informing the user about the min Teff required to        
c         neglect occultation by the central source.                                
c     =======================================================================       
      IMPLICIT NONE                                                             
      CHARACTER*10 Tstrg                                                        
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
      INTEGER iL                                                                
      DOUBLE PRECISION qaux(npL), qaux2(npL), res1, res2, Te_min,               
     &        mx, Psitn, Planck, xP                                             

c     Estimate min Teff required to neglect occultation (eq.(5) in Manual):     
      DO iL = 1, nL                                                             
        qaux(iL) = SigmaA(1,iL) * Us(iL,1) / lambda(iL)                         
        xP = 14400.0 / Tsub(1) / lambda(iL)                                     
        qaux2(iL) = SigmaA(1,iL) * Planck(xP) / lambda (iL)                     
      END DO                                                                    
      CALL Simpson(npL,1,nL,lambda,qaux,res1)                                   
      CALL Simpson(npL,1,nL,lambda,qaux2,res2)                                  
c     Approximate Psi for opt.thin case:                                        
      Psitn = res1/res2                                                         
      mx = Tsub(1)*sqrt(sqrt(4./Psitn))                                         
      IF(Tsub(1).LT.mx) THEN                                                    
       Te_min = 2. * mx                                                         
      ELSE                                                                      
       Te_min = 2. * Tsub(1)                                                    
      END IF                                                                    
      CALL GetFS(Te_min,0,1,Tstrg)                                              
      write(12,*)                                                               
     & ' ====================================================  '                
      write(12,*)                                                               
     & ' For compliance with the point-source assumption, the'                  
      write(12,*)                                                               
     & ' following results should only be applied to sources '                  
      write(12,'(a37,a5,a3)')                                                   
     & '  whose effective temperature exceeds ',Tstrg, ' K.'                    
      write(12,*)                                                               
     & ' ===================================================='                  

      RETURN                                                                    
      END                                                                       


      SUBROUTINE Pgrid(pstar,iPstar,error)                                      
c     =======================================================================       
c         After having the Y grid, generate the P grid (impact parameters)          
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
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER i, iP, j, error, Naux, iPstar, NinsLoc                            
      DOUBLE PRECISION pstar, delP                                              

      error = 0                                                                 
      P(1) = 0.0                                                                
      Naux = 1                                                                  
      iYfirst(Naux) = 1                                                         
      YPequal(Naux) = 1                                                         
c     grid within the cavity (points concentrated towards p=1, SQRT!)           
      nPcav = Ncav + Naux                                                       
      DO iP = Naux+1, nPcav                                                     
        P(iP) = P(Naux) + (1.0-P(Naux))*sqrt((iP-Naux)/(Ncav+1.0))              
        iYfirst(iP) = 1                                                         
        YPequal(iP) = 1                                                         
      END DO                                                                    
c     insert two rays if the stellar disk size is finite                        
      IF (pstar.GT.0.0) THEN                                                    
        CALL INSERT(pstar,iP,iPstar)                                            
        DO i = nPcav+1, iP                                                      
          iYfirst(i) = 1                                                        
          YPequal(i) = 1                                                        
        END DO                                                                  
        nPcav = iP                                                              
        ELSE                                                                    
        iP = nPcav                                                              
        iPstar = 0                                                              
      ENDIF                                                                     
c     impact parameters for the envelope (tangential to the radial grid)        
      DO i =1, nY-1                                                             
        Plast(i) = iP + 1                                                       
c       points close to the inner cavity should have more impact param.         
        NinsLoc = Nins                                                          
c       if tau is large add a few more points to the impact param. grid         
        IF (TAUmax.GT.10) THEN                                                  
          IF (i.LE.5) NinsLoc = Nins + 1                                        
        END IF                                                                  
        IF (TAUmax.GT.100) THEN                                                 
          IF (i.LE.10) NinsLoc = Nins + 1                                       
          IF (i.LE.5) NinsLoc = Nins + 2                                        
        END IF                                                                  
c       increment in P                                                          
        delP = (Y(i+1) - Y(i)) / (NinsLoc+1.0)                                  
        DO j = 1, NinsLoc+1                                                     
          iP = iP + 1                                                           
          iYfirst(iP) = i                                                       
          IF (j.eq.1) THEN                                                      
            YPequal(iP) = 1                                                     
            ELSE                                                                
            YPequal(iP) = 0                                                     
          END IF                                                                
          P(iP) = Y(i) + (j-1)*delP                                             
        END DO                                                                  
      END DO                                                                    
c     number of points for the P grid                                           
      nP = iP + 1                                                               
      P(nP) = Y(nY)                                                             
      iYfirst(nP) = nY                                                          
      YPequal(nP) = 1                                                           
      Plast(nY) = nP                                                            
c     check that the P grid is not too large                                    
      IF (nP.GT.npP) THEN                                                       
        error = 2                                                               
        CALL LINE(0,2,12)                                                       
        write(12,'(a)')' This model terminated because needed accuracy'         
        write(12,'(a,i6)')' results in too many points in P grid:',             
     &                    nP                                                    
        write(12,'(a,i4,a)')'          (See PARAMET.inc, npP =',npP,')'         
        CALL LINE(0,2,12)                                                       
        iERROR = iERROR + 1                                                     
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      DOUBLE PRECISION FUNCTION Planck(x)                                       
c     =======================================================================       
c     This function evaluates the Planck function multiplied by wavelength          
c     and normalized by sigma*T^4/Pi.                      [Z.I., Mar. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      DOUBLE PRECISION x                                                        

      IF (x.GT.100) THEN                                                        
        Planck = 0.0                                                            
        ELSE                                                                    
        IF (x.LT.0.00001) THEN                                                  
          Planck = 0.154 * x**3.0                                               
          ELSE                                                                  
          Planck = 0.154 * x**4.0 / (dexp(x) - 1)                               
        END IF                                                                  
      END IF                                                                    

      RETURN                                                                    
      END                                                                       

                                                                                

      SUBROUTINE RADTRANSF(pstar,iPstar,TAUlim,nG,ETAzp,FbolOK,deviat,          
     &                     error,iterFbol,model)                                
c     =======================================================================       
c     This subroutine solves the continuum radiative transfer problem for a         
c     spherically symmetric envelope.                      [Z.I., Nov. 1995]        
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
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER error, iPstar, nG, Conv, iter, FbolOK, iterFbol, iaux, iY,        
     &        Fconv, Uconv, BolConv, model, itnum, itlim                        
      DOUBLE PRECISION pstar, ETAzp(npP,npY),mat0(npL,npY,npY), maxFerr,        
     &       mat1(npL,npY,npY), Em(npL,npY),  alpha(npG,npY), deviat,           
     &       Uold(npL,npY), fbolold(npY), mifront(npL,npP,npY), dmaxU,          
     &       miback(npL,npP,npY), fsbol(npY), fdbol(npY), TAUlim, dmaxF         

c     generate, or improve, or do not touch the Y and P grids                   
      IF (iterETA.EQ.1.OR.iterFbol.GT.1) THEN                                   
        IF (iterFbol.EQ.1) THEN                                                 
c         first time generate grids                                             
          CALL SetGrids(pstar,iPstar,error,TAUlim)                              
          IF (error.NE.0) goto 999                                              
          IF (iVerb.EQ.2) write(*,*) 'Done with SetGrids'                       
        ELSE                                                                    
c         or improve the grid from previous iteration                           
          CALL ChkFlux(fBol,accuracy,iaux,error,ETAzp)                          
          IF (error.NE.0) goto 999                                              
c         generate new impact parameter grid                                    
c         increase the number of rays through the cavity                        
          IF (Ncav.LT.80) THEN                                                  
            Ncav = 2 * Ncav                                                     
            IF (iX.NE.0) write(18,'(a20,i3)')' Ncav increased to:',Ncav         
          END IF                                                                
c         increase the number of rays per y-grid interval                       
          IF (iterFbol.EQ.3.AND.Nins.EQ.2) THEN                                 
            Nins = Nins + 1                                                     
            IF (iX.NE.0) write(18,'(a20,i3)')' Nins increased to:',Nins         
          END IF                                                                
          CALL Pgrid(pstar,iPstar,error)                                        
c         if P grid is not OK end this model                                    
          IF (error.NE.0) goto 999                                              
          IF (iX.NE.0) THEN                                                     
            write(18,'(a23,i3)')' Y grid improved, nY =',nY                     
            write(18,'(a23,i3)')'                  nP =',nP                     
            write(18,'(a23,i3)')'                Nins =',Nins                   
            write(18,'(a23,i3)')'                Ncav =',Ncav                   
          END IF                                                                
        END IF                                                                  
      ELSE                                                                      
        IF (iX.NE.0) write(18,*)' Using same Y and P grids'                     
      END IF                                                                    
c     generate spline coefficients for ETA                                      
      CALL setupETA                                                             
c     evaluate ETAzp                                                            
      CALL getETAzp(ETAzp)                                                      
c     generate albedo through the envelope                                      
      CALL getOmega(nG)                                                         
c     generate stellar moments                                                  
      CALL Star(pstar,ETAzp,error)                                              
      IF (iVerb.EQ.2) write(*,*) 'Done with Star'                               
c     issue a message in fname.out about the condition for neglecting           
c     occultation:                                                              
      IF(model.eq.1.AND.iterFbol.eq.1.AND.iterEta.eq.1) CALL OccltMSG           
c     generate the first approximations for Td and alpha                        
      CALL InitTemp(ETAzp,nG,alpha)                                             
      IF (iVerb.EQ.2) write(*,*) 'Done with InitTemp'                           
c     find radiative transfer matrices                                          
      IF (iX.NE.0) write(18,*)' Calculating weight matrices'                    
      IF (iVerb.EQ.2) write(*,*) 'Calculating weight matrices'                  
      CALL Matrix(ETAzp,pstar,iPstar,mat0,mat1,mifront,miback)                  
      Conv = 0                                                                  
      iter = 0                                                                  
c     itlim is an upper limit on number iterations                              
      itlim = 10000                                                             
      IF (iX.NE.0) write(18,*)' Weight matrices OK, calculating Tdust'          
      IF (iVerb.EQ.2) write(*,*)' Weight matrices OK, calculating Tdust'        
      IF (iInn.eq.1) THEN                                                       
        write(38,'(a8,i5)') '    nY= ',nY                                       
        write(38,*) '    iter   maxFerr     dmaxU       dmaxF'                  
      END IF                                                                    
c     === Iterations over dust temperature =========                            
      DO WHILE (Conv.EQ.0.AND.iter.LE.itlim)                                    
        iter = iter + 1                                                         
c       find emission term                                                      
        CALL Emission(0,nG,Td,alpha,abund,Us,Em)                                
c       solve for Utot                                                          
        CALL Invert(nY,nL,mat0,Us,Em,omega,Utot,Uold,error)                     
        IF(error.NE.0) goto 999                                                 
c       find new Td and alpha, and check convergence                            
        CALL FindTemp(1,Utot,nG,Td,alpha)                                       
c       --------------------------------------                                  
c       every itnum-th iteration check convergence:                             
          IF (iter.GT.80) THEN                                                  
             itnum = 10                                                         
            ELSE                                                                
             itnum = 6                                                          
          END IF                                                                
c       first find 'old' flux (i.e. in the previous iteration)                  
        IF (MOD(iter+1,itnum).EQ.0) THEN                                        
          CALL Multiply(1,npY,nY,npL,nL,mat1,Utot,omega,0,fs,fds)               
          CALL Multiply(0,npY,nY,npL,nL,mat1,Em,omega,0,fs,fde)                 
          CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                               
c         find bolometric flux                                                  
          CALL Bolom(ftot,fbolold)                                              
        END IF                                                                  
        IF (MOD(iter,itnum).EQ.0) THEN                                          
c         first calculate total flux                                            
          CALL Multiply(1,npY,nY,npL,nL,mat1,Utot,omega,0,fs,fds)               
          CALL Multiply(0,npY,nY,npL,nL,mat1,Em,omega,0,fs,fde)                 
          CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                               
c         find bolometric flux                                                  
          CALL Bolom(ftot,fbol)                                                 
c         check convergence of bolometric flux                                  
          CALL Converg1(nY,accFbol,dynrange,fbolold,fbol,Fconv,dmaxF)           
c         check convergence of energy density                                   
          CALL Converg2(nY,nL,accConv,dynrange,Uold,Utot,Uconv,dmaxU)           
c         find maximal fbol error                                               
          CALL FindErr(fbol,maxFerr,nY)                                         
c         ------  printout of errors and convergence with iter.(inner flag): -------    
          IF(iInn.EQ.1) THEN                                                    
            write(38,'(i7,1p,3e12.4)') iter,maxFerr,dmaxU,dmaxF                 
          END IF                                                                
c         --------------------------------------------------------------                
          IF (maxFerr.LE.accuracy) THEN                                         
            BolConv = 1                                                         
          ELSE                                                                  
            BolConv = 0                                                         
          END IF                                                                
c         total criterion for convergence: Utot must converge, and ftot         
c         must either converge or have the required accuracy                    
          IF (Uconv*(Fconv+BolConv).GT.0) Conv = 1                              
        END IF                                                                  
c       --------------------------------------                                  
      END DO                                                                    
c     === The End of Iterations over Td ===                                     
      IF (iX.NE.0) THEN                                                         
        IF (iter.LT.itlim) write(18,*)                                          
     &     ' Convergence achieved, number of'                                   
        write(18,'(a34,i4)')                                                    
     &     ' iterations over energy density: ',iter                             
        write(18,'(a30,1p,e8.1)')                                               
     &     ' Flux conservation OK within:',maxFerr                              
        IF (iter.GE.itlim) THEN                                                 
          CALL MSG(2)                                                           
          iWARNING = iWARNING + 1                                               
        END IF                                                                  
      END IF                                                                    
c     calculate the emission term for the converged Td                          
      CALL Emission(0,nG,Td,alpha,abund,Us,Em)                                  
c     calculate flux                                                            
      CALL Multiply(1,npY,nY,npL,nL,mat1,Utot,omega,0,fs,fds)                   
      CALL Multiply(0,npY,nY,npL,nL,mat1,Em,omega,0,fs,fde)                     
      CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                                   
      CALL Bolom(ftot,fbol)                                                     
c     check whether, and how well, is bolometric flux conserved                 
      CALL ChkBolom(fbol,accuracy,deviat,FbolOK)                                
c     calculate additional output quantities                                    
c     1) energy densities                                                       
      CALL Multiply(1,npY,nY,npL,nL,mat0,Utot,omega,0,Us,Uds)                   
      CALL Multiply(0,npY,nY,npL,nL,mat0,Em,omega,0,Us,Ude)                     
      CALL Add(npY,nY,npL,nL,Us,Uds,Ude,Uchck)                                  
      CALL Bolom(Utot,Ubol)                                                     
      CALL Bolom(Uchck,UbolChck)                                                
c     2) scaled radial optical depth, tr                                        
      DO iY = 1, nY                                                             
        tr(iY) = ETAzp(1,iY) / ETAzp(1,nY)                                      
      END DO                                                                    
c     3) calculate intensity (at the outer edge) if required                    
      IF(iC.NE.0) THEN                                                          
       CALL FindInt(nG,alpha,ETAzp)                                             
       IF (iVerb.EQ.2) write(*,*) 'Done with FindInt'                           
      END IF                                                                    
c     if needed convolve intensity with the PSF                                 
      IF (iPSF.NE.0) THEN                                                       
        CALL Convolve(IntOut)                                                   
        IF (iVerb.EQ.2) write(*,*) 'Done with Convolve'                         
      END IF                                                                    
c     if needed find the visibility function                                    
      IF (iV.NE.0) THEN                                                         
        CALL Visibili(IntOut)                                                   
        IF (iVerb.EQ.2) write(*,*) 'Done with Visibili'                         
      END IF                                                                    
c     ============ if the inner flag iInn=1:  =========                             
      IF(iX.GE.1 .AND. iInn.EQ.1) THEN                                          
        CALL Bolom(fs,fsbol)                                                    
        CALL ADD2(fds,fde,fdbol,nY)                                             
        write(18,'(a11,1p,E11.3)')'   TAUfid =',TAUfid                          
        write(18,'(a11,1p,E11.3)')'  MaxFerr =',maxFerr                         
        write(18,*)                                                             
     &'     tr      fbol       fsbol      fdbol       Ubol  '                   
       DO iY = 1, nY                                                            
        write(18,'(1p,5E11.3)') tr(iY), fbol(iY), fsbol(iY),                    
     &                          fdbol(iY), Ubol(iY)                             
       END DO                                                                   
      END IF                                                                    
c     =====================                                                         

999   RETURN                                                                    
      END                                                                       


      SUBROUTINE SetGrids(pstar,iPstar,error,TAU)                               
c     =======================================================================       
c     Sets the Y and P grids based on GrayBody flux conservation.                   
c                                                         [MN & ZI, July'96]        
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
      INTEGER error, iPstar, dummy                                              
      DOUBLE PRECISION pstar, Ugb(npY), fgb(npY), ETAzp(npP,npY),albedo,        
     &       aux, faccs, TAU, accur, delTAUin                                   

c     store the default value for delTAUsc and facc                             
      faccs = facc                                                              
      delTAUin = delTAUsc                                                       
c     change the delTAUsc seed for the initial Y grid if TAU is large           
      IF (TAU.LT.1.0) delTAUsc = delTAUin * 2.0                                 
      IF (TAU.GE.1.0.and.TAU.LT.5.0) delTAUsc = delTAUin * 1.5                  
      IF (TAU.EQ.5.0) delTAUsc = delTAUin                                       
      IF (TAU.GT.5.0.and.TAU.LT.10.0) delTAUsc = delTAUin / 1.2                 
      IF (TAU.GE.10.0.and.TAU.LT.20.0) delTAUsc = delTAUin / 1.3                
c     The grid is set with TAU=min{TAUlim,TAUmax}, so the lines below are obsol
c     IF (TAU.GE.20.0.and.TAU.LT.30.0) delTAUsc = delTAUin / 1.4               
c     IF (TAU.GE.30.0.and.TAU.LT.50.0) delTAUsc = delTAUin / 1.5               
c     IF (TAU.GE.50.0) delTAUsc = delTAUin / 2.0                               
c     for steep density distributions (RDW, including analyt.approximation):    
      IF (denstyp.EQ.4.OR.RDW) delTAUsc = delTAUin / 1.2               
c     change the facc seed for the initial Y grid if Yout is very small         
      IF (Yout.LT.1000.0) facc = dsqrt(faccs)                                   
      IF (Yout.LT.100.0) facc = dsqrt(facc)                                     
      IF (Yout.LT.10.0) facc = dsqrt(facc)                                      
      IF (Yout.LT.2.0) facc = dsqrt(facc)                                       
      IF (Yout.LT.1.2) facc = dsqrt(facc)                                       
      IF (Yout.LT.1.05) facc = dsqrt(facc)                                      
                                                                                
      albedo = 1.0                                                              
      aux = 1.0                                                                 
c     generate initial grids                                                    
      CALL Ygrid(pstar,iPstar,error)                                            
      IF (error.NE.0) goto 101                                                  
      CALL Pgrid(pstar,iPstar,error)                                            
      IF (error.NE.0) goto 101                                                  
      IF (iX.GE.1) THEN                                                         
         write(18,'(a24,i3)')' Y grid generated, nY =',nY                       
         write(18,'(a24,i3)')'                   nP =',nP                       
         write(18,'(a24,i3)')'                 Nins =',Nins                     
         write(18,'(a24,i3)')'                 Ncav =',Ncav                     
      END IF                                                                    
c     solve for gray body (i.e. pure scattering)                                
      CALL GrayBody(albedo,TAU,Ugb,fgb)                                         
      IF (iVerb.EQ.2) write(*,*) 'Done with GrayBody'                           
c     find the max deviation of fgb (FindRMS called with flag 1)                
      CALL FindRMS(1,fgb,aux,accur,nY)                                          
      IF (iX.GE.1) THEN                                                         
        IF (accur.GT.accuracy) THEN                                             
          write(18,'(a25)')' Grids need improvement:'                           
          write(18,'(a29,1p,e10.3)')                                            
     &                          '                   fTot(nY):',fgb(nY)          
          write(18,'(a29,1p,e10.3)')'      Single wavelength TAU:',TAU          
          write(18,'(a29,1p,e10.3)')                                            
     &                         '          Required accuracy:',accuracy          
        END IF                                                                  
        write(18,'(a29,1p,e10.3)')' Single wavelength accuracy:',accur          
      END IF                                                                    
      IF(accur.GT.accuracy) THEN                                                
c       ChkFlux checks the bolometric flux conservation for the given           
c       grid and decreases the step if conservation is not satisfactory         
        dummy = 5                                                               
        CALL ChkFlux(fgb,accuracy,dummy,error,ETAzp)                            
        IF (error.NE.0) goto 101                                                
c       dummy=5 means everything was fine in ChkFlux                            
        IF (dummy.EQ.5) THEN                                                    
          IF (iX.GE.1) write(18,'(a23,i3)')' Y grid improved, nY =',nY          
c         generate new impact parameter grid                                    
          CALL Pgrid(pstar,iPstar,error)                                        
c         if P grid is not OK end this model                                    
          IF (error.NE.0) goto 101                                              
        ELSE                                                                    
          IF (iX.GE.1) THEN                                                     
            write(18,'(a59,i3)')                                                
     &      ' Although single wavelength accuracy was not satisfactory,'        
            write(18,'(a56,i3)')                                                
     &      ' Y grid could not be improved because npY is too small.'           
            write(18,'(a58,i3)')                                                
     &      ' Continuing calculation with a hope that it will be fine.'         
          END IF                                                                
        END IF                                                                  
      END IF                                                                    
c     return the default value for facc                                         
101    facc = faccs                                                             
       delTAUsc = delTAUin                                                      

      RETURN                                                                    
      END                                                                       

      SUBROUTINE solve(model,nG,error,ETAzp)                                    
c     =======================================================================       
c     This subroutine solves the continuum radiative transfer problem for a         
c     spherically symmetric envelope.                      [Z.I., Nov. 1995]        
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
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER model,error, iPstar, iterFbol, nG, FbolOK, EtaOK, iY              
      DOUBLE PRECISION pstar, ETAzp(npP,npY), TAUlim, deviat                    

      IF(iInn.eq.1) THEN                                                        
        write(38,*)'============================================='              
        write(38,'(a7,i5)') ' model= ',model                                    
        write(38,*)'============================================='              
      END IF                                                                    
      IF (iX.NE.0) THEN                                                         
       CALL LINE(0,2,18)                                                        
       write(18,'(a7,i3,a20)') ' model ',model,'  RUN-TIME MESSAGES '           
       CALL LINE(0,1,18)                                                        
       write(18,*)'===========  STARTING SOLVE  ==========='                    
      END IF                                                                    
      error = 0                                                                 
      IF(denstyp.NE.0) THEN                                                     
c      Solve for spherical envelope:                                             
c      temporarily the star is approximated by a point source                   
       pstar = 0.0                                                              
       iPstar = 1                                                               
c      select optical depth for the grid calculation                            
c      based on dynamical range                                                 
       TAUlim = 0.5*log(1.0/dynrange)                                           
c      if actual maximal tau is smaller use that value                          
       IF (TAUmax.LT.TAUlim) TAUlim = TAUmax                                    
c      counter over ETA (for radiatively driven winds only)                     
       iterETA = 0                                                              
       EtaOK = 0                                                                
c      iterations over ETA                                                      
       DO WHILE (EtaOK.EQ.0)                                                    
         iterETA = iterETA + 1                                                  
         IF (iX.NE.0.AND.RDW) THEN                                              
           write(18,*)'----------------------------------------'                
           write(18,*)' ',iterETA,'. iteration over ETA'                        
         END IF                                                                 
         IF (iVerb.EQ.2.AND.RDW)                                                
     &          write(*,*) ' ',iterETA,'. iteration over ETA'                   
c        counter for iterations over bolometric flux conservation               
         iterFbol = 0                                                           
         FbolOK = 0                                                             
         DO WHILE (FbolOK.EQ.0)                                                 
           iterFbol = iterFbol + 1                                              
           IF (iX.NE.0) THEN                                                    
            write(18,*)'  ',iterFbol,'. iteration over Fbol'                    
           END IF                                                               
           IF (iVerb.EQ.2) write(*,*) iterFbol,'. iteration over Fbol'          
c          solve the radiative transfer problem                                 
           CALL RADTRANSF(pstar,iPstar,TAUlim,nG,ETAzp,FbolOK,deviat,           
     &                   error,iterFbol,model)                                  
           IF (iVerb.EQ.2) write(*,*) 'Done with RadTransf'                     
c          error.EQ.3 : file with stellar spectrum not available                
           IF (error.EQ.3) goto 999                                             
c          error.EQ.5 : Singular matrix                                         
           IF (error.EQ.5) goto 999                                             
c          error.EQ.6 : Eta exceeds limitations                                 
           IF (error.EQ.6) goto 999                                             
c          error.EQ.2 : P grid was not produced                                 
           IF (error.EQ.2.AND.iterFbol.EQ.1.AND.iterETA.EQ.1) THEN              
c           if this is the first calculation end this model                     
            iERROR = iERROR + 1                                                 
            goto 999                                                            
           ELSE                                                                 
c           if this is a higher iteration use previous solution                 
            IF (error.EQ.2) THEN                                                
              IF (iX.NE.0.AND.iterFbol.GT.1) THEN                               
              write(18,*)' ======= IMPORTANT WARNING ======== '                 
              write(18,*)' In trying to conserve Fbol reached'                  
              write(18,*)' the limit for grid sizes. Flux is '                  
              write(18,'(a,1p,e9.3)')'  conserved to within ', deviat           
              write(18,*)' Treat all results with caution!'                     
              END IF                                                            
              IF (iX.NE.0.AND.iterFbol.EQ.1) THEN                               
              write(18,*)' ======== IMPORTANT  WARNING ======== '               
              write(18,*)' In trying to converge on ETA reached'                
              write(18,*)' the limit for grid sizes. Flux is '                  
              write(18,'(a,1p,e9.3)')'  conserved to within ', deviat           
              write(18,*)' Treat all results with caution!'                     
              END IF                                                            
              error = 0                                                         
              FbolOK = 2                                                        
              iWARNING = iWARNING + 1                                           
            END IF                                                              
           END IF                                                               
c          just in case...                                                      
           IF (error.EQ.1) THEN                                                 
            IF (iX.NE.0) THEN                                                   
              write(18,*)' *********  FATAL ERROR  *********'                   
              write(18,*)' * Something was seriously wrong *'                   
              write(18,*)' * Contact Z. Ivezic, M. Elitzur *'                   
              write(18,*)' *********************************'                   
            END IF                                                              
            iERROR = iERROR + 1                                                 
            goto 999                                                            
           END IF                                                               
c          if Fbol not conserved try again with a finer grid                    
           IF (FbolOK.EQ.0.AND.iterFbol.LT.10.AND.iX.NE.0) THEN                 
            write(18,*)'  ******** MESSAGE from SOLVE ********'                 
            write(18,*)'  Full solution does not conserve Fbol'                 
            write(18,*)'       Y       TAU/TAUtot        fbol'                  
            DO iY =1, nY                                                        
              write(18,'(1p,3e13.4)')Y(iY),                                     
     &                  ETAzp(1,iY)/ETAzp(1,nY),fbol(iY)                        
            END DO                                                              
            write(18,*)'  Trying again with finer grids'                        
           END IF                                                               
c          if could not conserve Fbol in 10 trials give it up                   
           IF (FbolOK.EQ.0.AND.iterFbol.GE.10) THEN                             
            IF (RDW) THEN                                                       
            IF (iX.NE.0) THEN                                                   
            write(18,*)' **********  WARNING from SOLVE  **********'            
            write(18,*)' Could not obtain required accuracy in Fbol'            
            write(18,'(a26,1p,e10.3)')'  The achieved accuracy is:',            
     &                                 deviat                                   
            write(18,*)' Will try to converge on the dynamics, but '            
            write(18,*)' treat all results with caution !!         '            
            write(18,*)' If accuracy<=0.01, or TAUmax>1000, this   '            
            write(18,*)' code probably cannot do it. Otherwise,    '            
            write(18,*)' please contact Z. Ivezic or M. Elitzur    '            
            write(18,*)' ******************************************'            
            END IF                                                              
            iWARNING = iWARNING + 1                                             
            FbolOK = 1                                                          
            ELSE                                                                
            IF (iX.NE.0) THEN                                                   
            write(18,*)' **********  WARNING from SOLVE  **********'            
            write(18,*)' Could not obtain required accuracy in Fbol'            
            write(18,'(a26,1p,e10.3)')'  The achieved accuracy is:',            
     &                                 deviat                                   
            write(18,*)' !!!!  Treat all results with caution  !!!!'            
            write(18,*)' If accuracy<=0.01, or TAUmax>1000, this   '            
            write(18,*)' code probably cannot do it. Otherwise,    '            
            write(18,*)' please contact Z. Ivezic or M. Elitzur    '            
            write(18,*)' ******************************************'            
            END IF                                                              
            iWARNING = iWARNING + 1                                             
            FbolOK = 2                                                          
            END IF                                                              
           END IF                                                               
c        end of loop over flux conservation                                     
         END DO                                                                 
c        for winds check if ETA has converged...                                
         IF ((RDW).AND.FbolOK.NE.2) THEN                                        
c         ptr(2) is specified in INPUT and controls converg. crit.              
          IF (ptr(2).LT.1.0e-6.AND.iterETA.GT.2)THEN                            
            EtaOK = 1                                                           
          ELSE                                                                  
            CALL WINDS(nG,EtaOK,ETAzp,ftot)                                     
          END IF                                                                
          IF (iterETA.GT.10.AND.EtaOK.EQ.0) THEN                                
            EtaOK = 2                                                           
            iWARNING = iWARNING + 1                                             
            IF (iX.NE.0) THEN                                                   
              write(18,*)' *********  WARNING  *********'                       
              write(18,*)' Could not converge on ETA in '                       
              write(18,*)' 10 iterations.'                                      
              write(18,*)' *********************************'                   
            END IF                                                              
          END IF                                                                
c         ...or otherwise finish right away                                     
         ELSE                                                                   
          EtaOK = 1                                                             
         END IF                                                                 
c      end of loop over ETA                                                     
       END DO                                                                   
      ELSE                                                                      
c      solve for slab case                                                      
       CALL SLBsolve(model,nG,error)                                            
c      error=4 means npY not large enough for oblique illumination grid         
       IF (error.eq.4) THEN                                                     
        CALL MSG(15)                                                            
        iWARNING = iWARNING + 1                                                 
        goto 999                                                                
       END IF                                                                   
      END IF                                                                    
c     analyze the solution and calculate some auxiliary quantities              
      CALL Analysis(model,ETAzp,error)                                          
      IF (iVerb.EQ.2) write(*,*) 'Done with Analysis'                           
      IF (iX.NE.0) THEN                                                         
        write(18,*)' ==== SOLVE successfully completed ====='                   
        write(18,*)' ======================================='                   
      END IF                                                                    

999   RETURN                                                                    
      END                                                                       


      SUBROUTINE STAR(pstar,ETAzp,error)                                        
c     =======================================================================       
c     This subroutine generates the stellar moments   [ZI, Nov'95; MN,Mar'99]       
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
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      CHARACTER*235 line                                                        
      INTEGER iY, iL, ios1, iLs, nLs, k, kstop, i, error, Nlambdam              
c     Nlambdam is the max number entries for a user supplied stellar spectrum   
      PARAMETER (Nlambdam = 10000)                                              
      DOUBLE PRECISION pstar, ETAzp(npP,npY), zeta, Bb, x, Planck, a, b,        
     &        lambdaS(Nlambdam), Llamstar(Nlambdam), Stellar(Nlambdam),         
     &        Lstar, fL(100), fpl(npL), llS(Nlambdam), lS(Nlambdam),            
     &        expow, dyn2, EMfunc                                               

c     for  startyp.GE.4 stellar spectrum is read from the file 'nameStar'       
c     it is unit=3 (1 is the input file)                                        
      IF ((startyp(1).GE.4).AND.(startyp(1).LE.6)) THEN                         
        open(3,ERR=998,file=nameStar(1),STATUS='OLD')                           
        rewind(3)                                                               
        read(3,'(a235)',ERR=998) line                                           
        read(3,'(a235)',ERR=998) line                                           
        read(3,'(a235)',ERR=998) line                                           
        ios1 = 0                                                                
        iLs = 0                                                                 
        DO WHILE (ios1.ge.0)                                                    
          read(3,*,END=900,ERR=998,iostat=ios1)a,b                              
          IF (ios1.ge.0) THEN                                                   
            iLs = iLs + 1                                                       
            lambdaS(iLs) = a                                                    
            IF (a.LE.0.0) goto 998                                              
c           it is assumed that Llamstar is L_Lambda, but...                     
c           if startyp.EQ.4 then file gives lambda*L_lambda                     
            IF (startyp(1).EQ.4) Llamstar(iLs) = b / a                          
c           if startyp.EQ.5 then file gives L_lambda                            
            IF (startyp(1).EQ.5) Llamstar(iLs) = b                              
c           if startyp.EQ.6 then file gives Lnu=lambda**2*L_lambda              
            IF (startyp(1).EQ.6) Llamstar(iLs) = b / a / a                      
           END IF                                                               
        END DO                                                                  
900     close(3)                                                                
        IF (iLs.LT.2) goto 998                                                  
        nLs = iLs                                                               
c       if input wavelengths in descending order turn them around               
        IF (lambdaS(1).GT.lambdaS(2)) THEN                                      
          DO iLs = 1, nLs                                                       
            llS(iLs) = lambdaS(iLs)                                             
            lS(iLs) = Llamstar(iLs)                                             
          END DO                                                                
          DO iLs = 1, nLs                                                       
            lambdaS(iLs) = llS(nLs+1-iLs)                                       
            Llamstar(iLs) = lS(nLs+1-iLs)                                       
          END DO                                                                
        END IF                                                                  
c       normalize stellar spectrum                                              
        CALL Simpson(Nlambdam,1,nLs,lambdaS,Llamstar,Lstar)                     
c       generate dimensionless stellar spectrum                                 
        DO iLs = 1, nLs                                                         
          Stellar(iLs) = lambdaS(iLs) * Llamstar(iLs) / Lstar                   
        END DO                                                                  
        ELSE                                                                    
c       if startyp.EQ.3 generate power-law spectrum                             
        IF (startyp(1).EQ.3) THEN                                               
          fL(1) = 1.0                                                           
          IF (Nlamtr(1).GT.1) THEN                                              
            DO i = 2, Nlamtr(1)                                                 
              fL(i) = fL(i-1) * (lamtr(1,i-1)/lamtr(1,i))**klam(1,i-1)          
            END DO                                                              
          END IF                                                                
          DO iL = 1, nL                                                         
            IF ((lambda(iL)-lamtr(1,1))*(lambda(iL)-                            
     &                              lamtr(1,Nlamtr(1)+1)).LE.0.0) THEN          
              kstop = 0                                                         
              k = 0                                                             
              DO WHILE (kstop.EQ.0)                                             
                k = k + 1                                                       
                IF (lambda(iL).GE.lamtr(1,k)) THEN                              
                  kstop = 1                                                     
                  fpl(iL) = fL(k) * (lamtr(1,k)/lambda(iL))**klam(1,k)          
                END IF                                                          
              END DO                                                            
            ELSE                                                                
              fpl(iL) = 0.0                                                     
            END IF                                                              
          END DO                                                                
        END IF                                                                  
      END IF                                                                    
      DO iY = 1, nY                                                             
c       effect of the star's finite size                                        
        IF (pstar.GT.0.0) THEN                                                  
          zeta = 2.*(1.-sqrt(1.-(pstar/Y(iY))**2.))*(Y(iY)/pstar)**2.           
        ELSE                                                                    
          zeta = 1.0                                                            
        END IF                                                                  
c       loop over wavelengths                                                   
        DO iL = 1, nL                                                           
          IF (startyp(1).EQ.1) THEN                                             
            Bb = 0.0                                                            
            DO k = 1, nBB(1)                                                    
              x = 14400.0 / lambda(iL) / Tbb(1,k)                               
              Bb = Bb + rellum(1,k)*Planck(x)                                   
            END DO                                                              
          ELSE IF (startyp(1).EQ.2) THEN                                        
              Bb = EMfunc(lambda(iL),Tbb(1,1),xSiO)                             
          ELSE IF (startyp(1).EQ.3) THEN                                        
              Bb = fpl(iL)                                                      
          ELSE IF (lambda(iL).GT.lambdaS(nLs)) THEN                             
c             for lambda longer than the longest entry in nameStar              
c             assume Rayleigh-Jeans tail                                        
              Bb = Stellar(nLs) * (lambdaS(nLs)/lambda(iL))**3.                 
          ELSE IF (lambda(iL).LT.lambdaS(1)) THEN                               
c             if shorter than the shortest assume 0                             
              Bb = 0.0                                                          
          ELSE                                                                  
c             if within limits interpolate                                      
              CALL LinInter(Nlambdam,nLs,lambdaS,Stellar,                       
     &                                             lambda(iL),iLs,Bb)           
          END IF                                                                
c   ----------  done with stellar spectrum ---------------                      
c         stellar part of en.density and flux                                   
          expow = ETAzp(1,iY)*TAUtot(iL)                                        
          Us(iL,iY) = Bb * zeta * dexp(-expow)                                  
c         flux                                                                  
          fs(iL,iY) = Bb * dexp(-expow)                                         
c         need to carry fsL to be compatible with the slab case                 
          fsL(iL,iY) = fs(iL,iY)                                                
          fsR(iL,iY) = 0.0                                                      
          dyn2 = dynrange*dynrange                                              
          IF (Us(iL,iY).LT.dyn2) Us(iL,iY) = 0.0                                
          IF (fs(iL,iY).LT.dyn2.AND.denstyp.NE.0) fs(iL,iY) = 0.0               
        END DO                                                                  
      END DO                                                                    
c     normalize stellar quantities with the stellar bolometric flux             
      CALL Bolom(fsL,fsbol)                                                     
      DO iY = 1, nY                                                             
        DO iL = 1, nL                                                           
           Us(iL,iY) = Us(iL,iY) / fsbol(1)                                     
           fs(iL,iY) = fs(iL,iY) / fsbol(1)                                     
           fsL(iL,iY) = fsL(iL,iY) / fsbol(1)                                   
        END DO                                                                  
      END DO                                                                    
      error = 0                                                                 
      goto 999                                                                  
998   write(12,*)' *** FATAL ERROR IN DUSTY! *************************'         
      write(12,*)' File with the spectral shape of external radiation:'         
      write(12,'(a2,a70)')'  ',nameStar(1)                                      
      write(12,*)' is missing or not properly formatted?!'                      
      write(12,*)' ***************************************************'         
      error = 3                                                                 

999   RETURN                                                                    
      END                                                                       


      SUBROUTINE trapzd2(a,b,s,n)                                               
c     =======================================================================       
c     This function integrates prescribed 8 functions from z=a to z=b with n        
c     divisions and stores the results to s(1..8). It is a heavily modified         
c     version of subroutine 'trapzd' (Num.Rec.'92).        [MN & ZI, Aug'96]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER it,iC,i,n,j                                                       
      DOUBLE PRECISION s(8),a,b,funcx(8),funca(8),funcb(8),del,sum(8),          
     &       tnm, x, ff, gp, gm                                                 

      IF (n.eq.1) then                                                          
c       calculate auxiliary functions at a and at b                             
        CALL TWOFUN(a,ff,gp,gm)                                                 
        funca(1) =  gm                                                          
        funca(5) =  gp                                                          
        DO iC = 2, 4                                                            
            funca(iC) = funca(iC-1) * ff                                        
            funca(4+iC) = funca(3+iC) * ff                                      
        END DO                                                                  
        CALL TWOFUN(b,ff,gp,gm)                                                 
        funcb(1) =  gm                                                          
        funcb(5) =  gp                                                          
        DO iC = 2, 4                                                            
          funcb(iC) = funcb(iC-1) * ff                                          
          funcb(4+iC) = funcb(3+iC) * ff                                        
        END DO                                                                  
c       calculate integrals for all 8 functions                                 
        DO i = 1, 8                                                             
          s(i) = 0.5*(b-a)*(funca(i)+funcb(i))                                  
        END DO                                                                  
      ELSE                                                                      
        it=2**(n-2)                                                             
        tnm=1.0*(it)                                                            
        del=(b-a)/tnm                                                           
        x=a+0.5*del                                                             
        DO i=1,8                                                                
          sum(i)=0.0                                                            
        END DO                                                                  
c       calculate contributions of all 'it' divisions                           
        DO j = 1, it                                                            
c         auxiliary functions at x                                              
          CALL TWOFUN(x,ff,gp,gm)                                               
c         generate (8) integrated functions at x                                
          funcx(1) = gm                                                         
          funcx(5) = gp                                                         
          DO iC = 2, 4                                                          
            funcx(iC) = funcx(iC-1) * ff                                        
            funcx(4+iC) = funcx(3+iC) * ff                                      
          END DO                                                                
          DO i=1,8                                                              
            sum(i)=sum(i)+funcx(i)                                              
          END DO                                                                
c         next x                                                                
          x=x+del                                                               
        END DO                                                                  
c       evaluate new value of the integral for all 8 cases                      
        DO i=1,8                                                                
          s(i)=0.5*(s(i)+(b-a)*sum(i)/tnm)                                      
        END DO                                                                  
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE TWOFUN(z,ff,gp,gm)                                             
c     =======================================================================       
c     This function evaluates auxiliary functions needed in trapzd2.                
c                                               [MN & ZI,Aug'96; MN,Sep'97]         
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
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
      INTEGER iLaux, iW1, iC                                                    
      DOUBLE PRECISION paux, w1,wL,etaloc, w, IntETA, pp, delTAUzp, z,          
     &        auxw, ff, gp, gm, gp1,gm1                                         
      COMMON /phi2/ paux, w1, wL, delTAUzp, iLaux, iW1                          

c     local radius                                                              
      w = dsqrt(paux*paux + z*z)                                                
      IF (w.LT.w1) w = w1                                                       
c     find local value for ETA function                                         
      etaloc = 0.0                                                              
      auxw = 1.                                                                 
      DO iC = 1, 4                                                              
        etaloc = etaloc + ETAcoef(iW1,iC)*auxw                                  
        auxw = auxw/w                                                           
      END DO                                                                    
c     ff, i.e. radial optical depth:                                            
      pp = 0.0                                                                  
      ff = IntETA(pp,iW1,wL,w)*TAUtot(iLaux)                                    
c     g functions:                                                              
      gp1 = dexp(IntETA(paux,iW1,w1,w)*TAUtot(iLaux)-delTAUzp)                  
      gm1 = dexp(-IntETA(paux,iW1,w1,w)*TAUtot(iLaux))                          
      gp = etaloc/w/w * gp1                                                     
      gm = etaloc/w/w * gm1                                                     

      RETURN                                                                    
      END                                                                       


      SUBROUTINE WEIGHTS(TAUaux,iP,iL,nZ,nY,alpha,beta,gamma,delta,             
     &                   wgp,wgm)                                               
c     =======================================================================       
c     This subroutine calculates weights wgp(iZ,iY) and wgm(iZ,iY) for              
c     integrations:                                                                 
c     INT(S(w)*exp(sign*ETAzp(iP,iZ')/w^2)dETAzp(iP,iZ')]                           
c     from ETAzp(iP,iZ) to ETAzp(iP,iZ+1), where w is local radius                  
c     corresponding to TAU(iP,iZ'), and sign=1 for wgp and -1 for wgm.              
c     Integrals are evaluated as:                                                   
c     INT = wg(iZ,1)*S(1) + wg(iZ,2)*S(2) + ... + wg(iZ,nY)*S(nY) with              
c     iZ=1..nZ-1. The method is based on approximation of S by cubic spline         
c     in radial optical depth given through matrices alpha, beta, gamma and         
c     delta (see MYSPLINE).                         [ZI,Dec'95;MN,Sep'97]           
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER iYfirst, YPequal, Plast                                           
      DIMENSION iYfirst(npP), YPequal(npP), Plast(npY)                          
      COMMON /Yfirst/ iYfirst, YPequal, Plast                                   
      INTEGER iP, iL, nZ, nY, iW, iZ, j                                         
      DOUBLE PRECISION alpha(npY,npY), beta(npY,npY), gamma(npY,npY),           
     &        delta(npY,npY), TAUaux(npL,npP,npY), K1p(npY),K2p(npY),           
     &        K3p(npY),K4p(npY), K1m(npY),K2m(npY),K3m(npY),K4m(npY),           
     &        wgp(npY,npY), wgm(npY,npY), waux                                  

c     generate integrals of 'TAUr**n'                                           
      CALL Kint4(TAUaux,iP,iL,nZ,K1p,K2p,K3p,K4p,K1m,K2m,K3m,K4m)               
c     loop over position on the line of sight                                   
      DO iZ = 1, nZ                                                             
        iW = iYfirst(iP) + iZ - 1                                               
c       loop over radial position                                               
        DO j = 1, nY                                                            
         IF (iZ.GT.1) THEN                                                      
          waux = alpha(iW-1,j)*K1p(iZ) + beta(iW-1,j)*K2p(iZ)                   
          wgp(iZ,j)=waux + gamma(iW-1,j)*K3p(iZ)+delta(iW-1,j)*K4p(iZ)          
         ELSE                                                                   
          wgp(1,j) = 0.0                                                        
         END IF                                                                 
         IF (iZ.LT.nZ) THEN                                                     
          wgm(iZ,j) = alpha(iW,j)*K1m(iZ) + beta(iW,j)*K2m(iZ)                  
          wgm(iZ,j) = wgm(iZ,j)+gamma(iW,j)*K3m(iZ)+delta(iW,j)*K4m(iZ)         
         ELSE                                                                   
          wgm(nZ,j) = 0.0                                                       
         END IF                                                                 
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE Ygrid(pstar,iPstar,error)                                      
c     =======================================================================       
c     This subroutine generates the radial (Y) grid.                                
c     Ncav is desired number of rays through the inner cavity (<=20), denstype      
c     is type of density law (see input file). Yout is the relative thickness,      
c     Yout=rout/r1. delTausc is maximum increase of the scaled optical depth        
c     tau/tauTot between two successive radial grid points and facc is the          
c     maximum ratio of their y coordinates. EtaRat is initialized in SUB Input      
c     and limits the ratio of two successive Eta(iY), in order to prevent           
c     nonphysical results in case of steep density distribution (pow>2). pstar is   
c     the impact parameter for star (0<=p<=1). First few points are prescribed      
c     arbitrarily. error is set to 1 if desired delTausc and facc result in too     
c     many points (> npY from userpar.inc). This subroutine calls function ETA      
c     which evaluates the normalized density profile.    [ZI, Nov'95; MN,Sep'99]    
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
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER i, istop,error,iPstar, iter, iYdummy, iY, itr, j, ir, jY,         
     &        irmax                                                             
      DOUBLE PRECISION pstar, dyR, dyT, ETA, aux, tausc, EtaTemp(npY),          
     &       ee, Yloc, nh, fac, dyF, yold, ynew, y1, delY1, Ymid,               
     &       TAUloc, rat, dY, Etamin, Etamax                                    

c     max number iter. over improving ratio of two Eta's                        
      irmax = 20                                                                
c     save old grid and values of Eta (important for denstyp = 5 or 6)          
      IF (nY.GT.0.AND.(RDW)) THEN                                               
        DO iY = 1, nY                                                           
          Yprev(iY) = Y(iY)                                                     
          EtaTemp(iY) = ETAdiscr(iY)                                            
        END DO                                                                  
        nYprev = nY                                                             
      END IF                                                                    
      pstar = pstar                                                             
      iPstar = iPstar                                                           
      y1 = 1.0                                                                  
      iter = 0                                                                  
101   error = 0                                                                 
      iter = iter + 1                                                           
c     resolve inner boundary in tau space                                       
c     from requiring TAU(2) = TAUmax*ETA(1)*(Y(2)-1) =~ 1                       
      Y(1) = y1                                                                 
      IF (TAUmax*ETA(y1).GT.2000.0) THEN                                        
        delY1 = 2.0 / TAUmax / ETA(y1)                                          
      ELSE                                                                      
        delY1 = 1.0 / TAUmax / ETA(y1)                                          
      END IF                                                                    
c     if very thin generate at least about 10-15 pts.                           
      IF (delY1.GT.0.01/ETA(y1)) delY1 = 0.01/ETA(y1)                           
c     do not push it over the limit of spline stability                         
      IF (delY1.LT.0.00005) delY1 = 0.00005                                     
      i = 1                                                                     
      istop = 0                                                                 
      DO WHILE (istop.NE.1)                                                     
        i = i + 1                                                               
        Y(i) = Y(i-1) + delY1                                                   
        IF (Y(i).GT.facc*Y(i-1)) THEN                                           
          Y(i) = facc * Y(i-1)                                                  
        ELSE                                                                    
          delY1 = 2. * delY1                                                    
        END IF                                                                  
        Ymid = dsqrt(Y(i)*Y(i-1))                                               
        TAUloc = TAUmax * (Y(i)-1.0) * ETA(Ymid)                                
        IF (Y(i).GE.1.01.OR.TAUloc.GT.10.0) istop = 1                           
c       in case of shell thinner than 1.01 comment the above line and           
c       uncomment the next line, unless it is a case of RDW at high TauV.       
c       These require more points near the origin and 1.01 is a better limit.   
c       IF (Y(i).GE.1.0001.OR.TAUloc.GT.10.0) istop = 1                         
      END DO                                                                    
      Ynew = Y(i)                                                               
c     some rule of thumb estimates for factor nh                                
      IF (TAUmax.GT.10000.0) THEN                                               
c       extreme taus                                                            
        nh = 15.0                                                               
        Ncav = 80                                                               
      ELSE IF (TAUmax.GT.2000.0) THEN                                           
c         huge taus                                                             
          nh = 10.0                                                             
          Ncav = 40                                                             
      ELSE IF (TAUmax.GT.400.0) THEN                                            
c         large taus                                                            
          nh = 8.0                                                              
          Ncav = 20                                                             
      ELSE                                                                      
c         normal' taus                                                          
          nh = 5.0                                                              
          Ncav = 10                                                             
      END IF                                                                    
c     very small taus                                                           
      IF (TAUmax.LT.10.0) nh = nh / 2.0                                         
      IF (TAUmax.LT.1.0) nh = nh / 2.0                                          
c     empirically: ~1/r needs more points for small y:                          
      IF ((pow-1.4)*(pow-0.6).LE.0.0) nh = nh * 1.5                             
c     same for for steep density distributions:                                 
      IF (denstyp.eq.4.OR.RDW) nh = nh * 1.5                           
      tausc = ETA(Y(i))+ETA(Y(i-1))                                             
      tausc = (Y(i)-Y(i-1)) * 0.5 * tausc                                       
      fac = dexp(-dlog(tausc) / nh)                                             
      istop = 0                                                                 
c     for broken power-laws                                                     
      IF (Ntr.GE.1) THEN                                                        
        itr = 1                                                                 
        DO j = 1, i                                                             
          IF (Y(j).GE.Ytr(itr)) itr = itr + 1                                   
        END DO                                                                  
      END IF                                                                    
c     generate the rest of Y grid points                                        
      DO WHILE (istop.NE.1)                                                     
        i = i + 1                                                               
        Yold = Ynew                                                             
c       find maximal increase in Y allowed by the ratio facc                    
        dyR = Yold * (facc-1.)                                                  
c       find maximal increase in Y allowed by delTausc                          
        dyT = delTausc / ETA(Yold)                                              
c       find maximal increase in Y allowed by the ratio of tausc                
        dyF = tausc*(fac-1.) / ETA(Yold)                                        
c       find new Y                                                              
c        Ynew = Yold + MIN(dyR,dyT,dyF)                                         
        dY = MIN(dyR,dyT,dyF)                                                   
c       Check if the max ratio btw. two Eta values is less than EtaRat          
c       and insert additional y-pts. where necessary. This prevents sharp       
c       drops in Utot(npL,npY) in case of steep Eta's [MN'99].                  
        DO ir = 1 , irmax                                                       
          Ynew = Yold + dY                                                      
          rat = ETA(Yold)/ETA(Ynew)                                             
          IF(rat.GE.1./EtaRat .AND. rat.LE.EtaRat) goto 10                      
          dY = 0.5*dY                                                           
        END DO                                                                  
        CALL MSG(16)                                                            
10      continue                                                                
        Y(i) = Ynew                                                             
c       make sure that all transition points are included in the grid           
        IF (Ntr.GE.1) THEN                                                      
         IF (Y(i).GE.Ytr(itr)) THEN                                             
           Y(i) = Ytr(itr)                                                      
           Ynew = Y(i)                                                          
           itr = itr + 1                                                        
         END IF                                                                 
        END IF                                                                  
        aux = ETA(Ynew)+ETA(Yold)                                               
        aux = (Ynew-Yold) * 0.5 * aux                                           
        tausc = tausc + aux                                                     
c       finish when Yout is reached                                             
        IF (Ynew.GE.Yout) istop = 1                                             
      END DO                                                                    
      Y(i) = Yout                                                               
      nY = i                                                                    
c     insert additional penultimate point to avoid spline oscillations          
      Y(nY+1) = Yout                                                            
      Y(nY) = dsqrt(Y(nY)*Y(nY-1))                                              
      nY = nY + 1                                                               
c     check that outer edge is well resolved in tau space                       
c     (important for flat ETAs)                                                 
      istop = 0                                                                 
      DO WHILE (istop.NE.1)                                                     
        IF ((Yout-Y(nY-1))*TAUmax*ETA(Yout).GT.1.0) THEN                        
          Y(nY+1) = Yout                                                        
          Y(nY) = dsqrt(Y(nY)*Y(nY-1))                                          
          nY = nY + 1                                                           
        ELSE                                                                    
          istop = 1                                                             
        END IF                                                                  
      END DO                                                                    
c     check dynamical range of Eta to avoid nonphysical results or code errors [
      Etamax = 0.                                                               
      Etamin = 1.e+20                                                           
      DO iY = 1, nY                                                             
       IF(ETA(Y(iY)).lt.Etamin) Etamin = ETA(Y(iY))                             
       IF(ETA(Y(iY)).gt.Etamax) Etamax = ETA(Y(iY))                             
       IF (ETA(Y(iY)).LT.1.e-12) THEN                                           
        IF (iX.GT.0) THEN                                                       
         write(18,*)'      Y          ETA  '                                    
         DO jY = 1, nY                                                          
           write(18,'(1p,2e12.3)') Y(jY),ETA(Y(jY))                             
         END DO                                                                 
        END IF                                                                  
        CALL MSG(17)                                                            
        error = 6                                                               
        iERROR = iERROR + 1                                                     
        goto 102                                                                
       END IF                                                                   
      END DO                                                                    
      IF ((Etamin/Etamax).LT.1.e-12) THEN                                       
       CALL MSG(18)                                                             
       error = 6                                                                
       iERROR = iERROR + 1                                                      
       goto 102                                                                 
      END IF                                                                    
c     check that the Y grid is not too large                                    
      IF (nY.GT.npY) THEN                                                       
        delTAUsc = delTAUsc * 1.5                                               
        iWARNING = iWARNING + 1                                                 
        IF (iX.NE.0) THEN                                                       
        IF (iter.EQ.1) call line(0,2,18)                                        
        write(18,'(a)')                                                         
     &  ' ****************   WARNING   *******************'                     
        write(18,'(a46,i3)')                                                    
     &             ' Initial delTAUsc resulted in too many points:',nY          
        write(18,'(a,i3,a)')' You need to increase npY in userpar.inc'          
        write(18,*)'Multiplying delTAUsc by 1.5 and trying again'               
        write(18,'(a14,1p,e10.3)')' New delTAUsc:',delTAUsc                     
        END IF                                                                  
        IF (iter.LT.5) THEN                                                     
          goto 101                                                              
        ELSE                                                                    
          IF (iX.NE.0) THEN                                                     
            write(18,'(a)')                                                     
     &            ' ****************  GIVING UP  *******************'           
            call line(0,2,18)                                                   
          END IF                                                                
          error = 2                                                             
          goto 102                                                              
        END IF                                                                  
      END IF                                                                    
c     intepolate ETAdiscr to new Y grid (for RDW (denstyp=5 or 6))              
c      write(18,*)' ***** From Ygrid *****'                                     
c      write(18,*)'      Y      ETAdiscr      '                                 
      DO iY = 1, nY                                                             
        Yloc = Y(iY)                                                            
        IF (iterETA.GT.1) THEN                                                  
          CALL LinInter(npY,nYprev,Yprev,EtaTemp,Yloc,iYdummy,ee)               
          ETAdiscr(iY) = ee                                                     
        ELSE                                                                    
          ETAdiscr(iY) = ETA(Yloc)                                              
        END IF                                                                  
c        write(18,'(1p,2e12.4)')Y(iY), ETAdiscr(iY)                             
      END DO                                                                    
c      write(18,*)' *************************'                                  

102   RETURN                                                                    
      END                                                                       