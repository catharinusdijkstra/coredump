c     =======================================================================       
c     These are the subroutines for a full dynamical calculation of             
c     radiatively driven winds.                              [MN, Mar'99]       
c     =======================================================================       
C     Table of Contents                                                         
C                                                                               
C     CALCETA                                                                   
C     DYNAMICS                                                                  
C     EMFUNC                                                                    
C     GAMMAFUN                                                                  
C     UFUN                                                                      
C     VRATFUN                                                                   
C     WINDS                                                                     
c     =======================================================================       
                                                                                

      SUBROUTINE CalcETA(Y,qF,u,vrat,Eta,tauF,nY)                               
c     =======================================================================       
c     Calculates the dimensionless density profile ETA(y) and tauF(y) when          
c     taking the drift into account                         [ZI & MN, Aug'96]       
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER iY, nY                                                            
      DOUBLE PRECISION qF(npY), Y(npY), u(npY), Eta(npY), tauF(npY),            
     &       F1(npY), F2(npY), vrat(npY), EtaINT, INT                           

c     generate ETA and its integral (normalization constant)                    
      DO iY = 1, nY                                                             
        F1(iY) = vrat(iY)/u(iY)/Y(iY)/Y(iY)                                     
      END DO                                                                    
      CALL SIMPSON(npY,1,nY,Y,F1,EtaINT)                                        
c     find tauF                                                                 
      DO iY = 1, nY                                                             
        Eta(iY) = F1(iY)/EtaINT                                                 
        F2(iY) = qF(iY)*ETA(iY)                                                 
        CALL SIMPSON(npY,1,iY,Y,F2,INT)                                         
        tauF(iY) = TAUfid*INT                                                   
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE DYNAMICS(nG,ETAzp,ftot)                                        
c     =======================================================================       
c     This subroutine finds the velocity structure of a radiatively driven          
c     wind.                                                                         
c     *** This version works for single size grains only ***                        
c                                                           [ZI & MN, Aug'96]       
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
      DOUBLE PRECISION ugas(npY), qF(npY), vrat(npG,npY), Gamma(npY),           
     &       I1, I2, I3, CMdot, Cve, CM, Cr1                                    
      COMMON /dyn/ ugas, qF, vrat, Gamma, I1, I2, I3, CMdot, Cve, CM,           
     &       Cr1                                                                
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER nG, iL, iY, iconv, itr, ETAconv, uconv                            
      DOUBLE PRECISION ETAold(npY), u(npY), uold(npY), resaux,tauF(npY),        
     &       Faux(npL), u1, eps, uacc, GammaMax, ETAzp(npP,npY),                
     &       ftot(npL,npY),Sigaux                                               

c     Temporary local rescaling (to reconcile with Zeljko's f-lae which         
c     were for <sigma>/<V>) [MN]                                                
      Sigaux = SigExfid/aveV                                                    
      IF (iX.GE.1) THEN                                                         
         write(18,*)' Doing Dynamics'                                           
      END IF                                                                    
      IF(iVerb.EQ.2) write(*,*)' Doing Dynamics'                                
c     so far it works for nG=1 only:                                            
      IF (nG.GT.1) THEN                                                         
        write(12,*)' **************************** '                             
        write(12,*)' Change dynamics sub to nG>1! '                             
        write(12,*)'       PROGRAM STOPPED        '                             
        write(12,*)' **************************** '                             
        stop                                                                    
      END IF                                                                    
c     accuracy for velocity convergence same as for Utot:                       
      uacc = accConv                                                            
c     find the flux averaged extinction efficiency, qF(y), and flux             
c     averaged optical depth, tauF(y) (with old ETA, i.e. it can be             
c     determined by using ETAzp, see below)                                     
      DO iY = 1, nY                                                             
c       generate auxiliary function for lambda integration:                     
        DO iL = 1, nL                                                           
          Faux(iL) = (SigmaA(1,iL)+SigmaS(1,iL))*ftot(iL,iY)/lambda(iL)         
        END DO                                                                  
        CALL Simpson(npL,1,nL,lambda,Faux,resaux)                               
c       qF scaled with Qext at lamfid                                           
        qF(iY) = resaux / SigExfid                                              
c       tauF                                                                    
        IF (iY.EQ.1) THEN                                                       
          tauF(iY) = 0.0                                                        
        ELSE                                                                    
          DO iL = 1, nL                                                         
           Faux(iL)=ETAzp(1,iY)*TAUtot(iL)*ftot(iL,iY)/lambda(iL)               
          END DO                                                                
          CALL Simpson(npL,1,nL,lambda,Faux,resaux)                             
          tauF(iY) = resaux                                                     
        END IF                                                                  
      END DO                                                                    
c     assign input parameters to local variables                                
      GammaMax = ptr(1)                                                         
      IF (denstyp.EQ.5) THEN                                                    
        eps = pow                                                               
        u1 =  tauF(nY)*eps*(1.-0.5*GammaMax)/(1.-eps)                           
        ELSE                                                                    
        u1 = pow                                                                
      END IF                                                                    
c     Initial approximation for u(y)                                            
      DO iY = 1, nY                                                             
        uold(iY) = u1 + tauF(iY)*(1.-0.5*GammaMax)                              
      END DO                                                                    
c     Initial approximation for v/vd                                            
      CALL vratFun(Sigaux,qF,uold,vrat,nY)                                      
c     Initial approximation for Gamma                                           
      CALL GammaFun(GammaMax,qF,uold,vrat,Y,Gamma,nY,I1,I2,I3)                  
      iconv = 0                                                                 
      itr = 0                                                                   
c     ITERATIONS until u and eta converge within uacc                           
      DO WHILE(iconv.NE.1)                                                      
        itr = itr + 1                                                           
c       Find the new value of u(y)                                              
        IF (denstyp.EQ.5)                                                       
     &      u1 = tauF(nY)*eps*(1.-Gamma(nY))/(1.-eps)                           
        CALL uFun(u1,tauF,Gamma,u,nY)                                           
c       Find the new ETAdiscr(y) and tauF(y)                                    
        CALL CalcETA(Y,qF,u,vrat,ETAdiscr,tauF,nY)                              
c       generate new v/vd                                                       
        CALL vratFun(Sigaux,qF,u,vrat,nY)                                       
c       generate new Gamma                                                      
        CALL GammaFun(GammaMax,qF,u,vrat,Y,Gamma,nY,I1,I2,I3)                   
c       check convergence of u and Eta                                          
        IF (itr.GT.1) THEN                                                      
          CALL ChkConv(nY,uacc,uold,u,uconv)                                    
          CALL ChkConv(nY,uacc,ETAold,ETAdiscr,ETAconv)                         
c         convergence required for both u(y) and ETA(y)                         
          iconv = ETAconv * uconv                                               
        END IF                                                                  
        IF (iconv.NE.1) THEN                                                    
          DO iY =1, nY                                                          
            uold(iY) = u(iY)                                                    
            ETAold(iY) = ETAdiscr(iY)                                           
          END DO                                                                
          IF (itr.GE.100) iconv = 1                                             
        ELSE                                                                    
          DO iY = 1, nY                                                         
            ugas(iY) = u(iY)                                                    
          END DO                                                                
        END IF                                                                  
      END DO                                                                    
      IF (iX.GE.1)                                                              
     &  write(18,'(a35,i3)')' Number of iterations to converge:',itr            

      RETURN                                                                    
      END                                                                       


      DOUBLE PRECISION FUNCTION EMfunc(lambda,Teff,xSiO)                        
c     This is modeled after subroutine engelke by M. Marengo. Here are his          
c     original comments:                                                            
c     =================================================================         
c     This subroutine computes a modified black body spectrum using an          
c     "Engelke" function (see Engelke 1992, AJ 104, 1248):                      
c     Bnu = Bnu(Tb) with Tb = 0.738*Teff*(1+79450/(lambda*Teff))**0.182         
c                                                                               
c     Molecular SiO absorption is modelled from the alpha Tau spectrum          
c     of Cohen et al. 1992, AJ 104,2030 with a 5th order polinomial,            
c     and added to the modified bb.                                             
c                                                                               
c     M. Marengo - mmarengo@cfa.harvard.edu - Sep 1998                          
c     =================================================================         
c                                                                               
c     This version makes use of the scaled quantities and Dusty's function          
c     Planck(x)                                                  [ZI, Feb 99]       
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER j                                                                 
      DOUBLE PRECISION lambda, Teff, xSiO, x, Planck, Teng, SiOc(6),            
     &                 lambda1, lambda2, SiO8m, SiOf                            

c     SiO fit data from Massimo:                                                   
c     Polinomial coeff for SiO absorption model (5th order),                    
c     wavelength interval in which to apply the absorption                      
c     and given absorption at 8 micron (to rescale for SiO)                     
      lambda1 =  7.8636                                                         
      lambda2 = 11.4280                                                         
      SiO8m = 1.0701447                                                         
      SiOc(1) = -300.43916                                                      
      SiOc(2) =  149.32134                                                      
      SiOc(3) =  -29.493280                                                     
      SiOc(4) =    2.9067144                                                    
      SiOc(5) =   -0.14304663                                                   
      SiOc(6) =    0.0028134070                                                 

c     Engelke's effective temperature                                           
      Teng = 0.738 * Teff * (1.0 + 79450.0/(lambda*Teff))**0.182                
      x = 14400.0 / lambda / Teng                                               
      EMfunc = (Teng/Teff)**4 * Planck(x)                                       
c     If lambda is in SiO region, compute and apply the SiO absorption          
      IF ((lambda-lambda1)*(lambda-lambda2).LT.0) THEN                          
          SiOf = 0.0                                                            
          DO j = 1, 6                                                           
             SiOf = SiOf + SiOc(j) * lambda**(1.0*j-1)                          
          END DO                                                                
          EMfunc = EMfunc / (1.0 + (SiOf-1)/(SiO8m-1)*xSiO*0.01)                
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE GammaFun(GammaMax,qF,u,vrat,Y,Gamma,nY,Int1,Int2,Int3)        
c     =======================================================================       
c     This subroutine finds the ratio of the gravitaional to the radiative          
c     force for the whole envelope Gamma^(-1). Note that this quantity is           
c     called Gamma in the code.                                  [MN, Aug'96]       
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER jY, nY                                                            
      DOUBLE PRECISION qF(npY), u(npY), vrat(npY), Gamma(npY), Y(npY),          
     &       Int2, Int3, Faux1(npY), Faux2(npY), Faux3(npY), GammaMax,          
     &       C, Int1                                                            

c     generate integrals                                                        
      DO jY = 1, nY                                                             
        Faux1(jY) = 1.0/u(jY)/Y(jY)/Y(jY)                                       
        Faux2(jY) = vrat(jY)*Faux1(jY)                                          
        Faux3(jY) = qF(jY)*Faux2(jY)                                            
        CALL SIMPSON(npY,1,jY,Y,Faux1,Int1)                                     
        CALL SIMPSON(npY,1,jY,Y,Faux2,Int2)                                     
        CALL SIMPSON(npY,1,jY,Y,Faux3,Int3)                                     
        IF (jY.GT.1) THEN                                                       
          Gamma(jY) = Int1 / Int3                                               
          ELSE                                                                  
          Gamma(jY) = 1.0 / vrat(1) / qF(1)                                     
        END IF                                                                  
      END DO                                                                    
c     find normalization constants                                              
      CALL FindMax(npY,1,nY,Gamma,C)                                            
      C = GammaMax / C                                                          
c     normalize Gamma                                                           
      DO jY = 1, nY                                                             
        Gamma(jY) = C * Gamma(jY)                                               
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE uFun(u1,tauF,Gamma,u,nY)                                       
c     =======================================================================       
c     Calculates the scaled gas velocity u(y) with given u1=u(1), tauF(y),          
c     and gravitational correction term Gamma^(-1) (named Gamma in the code).       
c                                                                [MN, Aug'96]       
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY,iY                                                             
      DOUBLE PRECISION u(npY), Gamma(npY), tauF(npY), u1                        

      DO iY = 1, nY                                                             
         u(iY) = u1 + tauF(iY) * (1.0-Gamma(iY))                                
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE vratFun(Qfid,qF,u,vrat,nY)                                     
c     =======================================================================       
c     Calculates the drift correction v(y)/vd(y).                [MN, Aug'96]       
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, iY                                                            
      DOUBLE PRECISION Qfid, qF(npY), u(npY), vrat(npG,npY)                     

      DO iY = 1, nY                                                             
         vrat(1,iY) = 1.0 / (1.0 + dsqrt(Qfid*qF(iY)/u(iY)))                    
      END DO                                                                    

      RETURN                                                                    
      END                                                                       
                                                                                

      SUBROUTINE WINDS(nG,EtaOK,ETAzp,ftot)                                     
c     =======================================================================       
c     This subroutine takes care of the interface between radiatively driven        
c     winds and radiative transfer.                              [ZI, Aug'95]       
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
      INTEGER nG, EtaOK, iY                                                     
      DOUBLE PRECISION ETAold(npY), ETAzp(npP,npY), ftot(npL,npY),              
     &       accETA                                                             

      DO iY = 1, nY                                                             
        ETAold(iY) = ETAdiscr(iY)                                               
      END DO                                                                    
c     find new ETA (ETAdiscr in density.inc) based on current flux ftot         
      CALL DYNAMICS(nG,ETAzp,ftot)                                              
c     check convergence (ptr(2) is specified in INPUT)                          
      accETA = ptr(2) * accuracy                                                
      CALL ChkConv(nY,accETA,ETAold,ETAdiscr,EtaOK)                             
      IF (iX.GE.1) THEN                                                         
        write(18,*)'     Y         ETAold      ETAnew      ratio'               
        DO iY = 1, nY                                                           
          accETA = ETAold(iY) / ETAdiscr(iY)                                    
          write(18,'(1p,4e12.5)')Y(iY),ETAold(iY),ETAdiscr(iY),accETA           
        END DO                                                                  
        IF (EtaOK.EQ.1) THEN                                                    
          write(18,*)' Convergence on Eta achieved'                             
        ELSE                                                                    
          write(18,*)' Convergence on Eta not achieved.'                        
          write(18,*)' Going to the next iteration.'                            
        END IF                                                                  
      END IF                                                                    
c     save Y to Yprev and nY to nYprev                                          
      DO iY = 1, nY                                                             
         Yprev(iY) = Y(iY)                                                      
      END DO                                                                    
      nYprev = nY                                                               

      RETURN                                                                    
      END                                                                       
