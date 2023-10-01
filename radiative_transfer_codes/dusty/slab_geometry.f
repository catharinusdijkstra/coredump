c     =========================================================================     
c              Below are the subroutines for calculation in slab geometry           
c                         arranged in alphabetical order                            
c     =========================================================================     
C         Table of Contents                                                         
C                                                                                   
C         SLBANALYT                                                                 
C         SLBDIFF                                                                   
C         EINT1                                                                     
C         EINT2                                                                     
C         EINT3                                                                     
C         EINT4                                                                     
C         EINT5                                                                     
C         SLBACC                                                                    
C         SLBGRAY                                                                   
C         SLBINIT                                                                   
C         SLBMAT                                                                    
C         SLBMISC                                                                   
C         SLBRADT                                                                   
C         SLBSOLVE                                                                  
C         SLBSTAR                                                                   
C         SLBTAU                                                                    
C         SLBTRACE                                                                  
C         SLBY                                                                      
c     =========================================================================     
                                                                                

      SUBROUTINE SLBAnalyt(TAU,L1,L2,L3,L4,N,NN)                                
c     =========================================================================     
c     Finds the integrals of the moments of the first exponential function E1   
c     Lk(Tau)=Int{t^(k-1) E1|t-Tau|}, k=1..4. The analytical expressions used   
c     are obtained integrating Lk by parts.                     [MN,Nov'98]     
c     =========================================================================     
      IMPLICIT NONE                                                             
      INTEGER i, j, N, NN                                                       
      DOUBLE PRECISION L1(NN,NN), L2(NN,NN), L3(NN,NN), L4(NN,NN), del,         
     &          TAU(NN), arg, arg1, Eint2, Eint3, Eint4, Eint5, Lt1,            
     &          Lt2, Lt3, Lt4, E2a1, E3a1, E4a1                                 

      DO i = 1, N                                                               
       DO j = 1, N                                                              
        L1(i,j) = 0.                                                            
        L2(i,j) = 0.                                                            
        L3(i,j) = 0.                                                            
        L4(i,j) = 0.                                                            
       END DO                                                                   
      END DO                                                                    
c     Calculate Lk(i)=Int_t^t1 {x**(k-1)E1|x-TAU(i)| dx}                        
      DO j = 1, N                                                               
        DO i = 1, N-1                                                           
          del = TAU(i+1)-TAU(i)                                                 
          arg1 = dabs(TAU(j)-TAU(i+1))                                          
          arg = dabs(TAU(j)-TAU(i))                                             
                                                                                
          E2a1 = Eint2(arg1)                                                    
          E3a1 = Eint3(arg1)                                                    
          E4a1 = Eint4(arg1)                                                    
          Lt1 = E2a1-Eint2(arg)                                                 
          Lt2 =-(E3a1-Eint3(arg))/del                                           
          Lt3 = E2a1 + 2./del/del*(E4a1-Eint4(arg))                             
          Lt4 =-3./del*E3a1 -6./del/del/del*                                    
     *                     (Eint5(arg1)-Eint5(arg))                             
          IF(i.LT.j) THEN                                                       
             L1(i,j) = 0.5*Lt1                                                  
             L2(i,j) = 0.5*(Lt2 + E2a1)                                         
             L3(i,j) = 0.5*(Lt3 - 2./del*E3a1)                                  
             L4(i,j) = 0.5*(Lt4 + E2a1 + 6./del/del*E4a1)                       
           ELSE                                                                 
            L1(i,j) = -0.5*Lt1                                                  
            L2(i,j) = 0.5*(Lt2 - E2a1)                                          
            L3(i,j) = 0.5*(-Lt3 - 2./del*E3a1 )                                 
            L4(i,j) = 0.5*(Lt4 - E2a1 - 6./del/del*E4a1)                        
          END IF                                                                
        END DO                                                                  
      END DO                                                                    
                                                                                
      DO i =1 , N                                                               
        DO j = 1, N                                                             
          IF(L1(i,j).LT.1.e-20) L1(i,j) = 0.                                    
          IF(L2(i,j).LT.1.e-20) L2(i,j) = 0.                                    
          IF(L3(i,j).LT.1.e-20) L3(i,j) = 0.                                    
          IF(L4(i,j).LT.1.e-20) L4(i,j) = 0.                                    
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
      SUBROUTINE SLBdifF(flag,om,grid,mat1,nL,nY,mat2,fp,fm)                    
c     =========================================================================     
c         Integration of U(lam,t)*E2|Tau-t| to get the diffuse flux.                
c         flag=1 is for scattered and flag=0 for emitted flux. [MN, Apr'98]         
c     =========================================================================     
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER iL, iY, j, nL, nY, flag                                           
      DOUBLE PRECISION mat1(npL,npY), mat2(npL,npY), grid(npL,npY),             
     &             om(npL,npY), fp(npL,npY), fm(npL,npY), TAU(npY),             
     &             faux(npY), Efact, fave, SUM, Eint3                           

      DO iL = 1, nL                                                             
        DO iY = 1, nY                                                           
          TAU(iY) = grid(iL,iY)                                                 
          IF (flag.EQ.1) THEN                                                   
             faux(iY) = om(iL,iY)*mat1(iL,iY)                                   
          ELSE                                                                  
             faux(iY) = (1.0-om(iL,iY))*mat1(iL,iY)                             
          END IF                                                                
        END DO                                                                  
                                                                                
c       Find f(+) (in arg Tau>t)                                                
        DO iY = 1, nY                                                           
          SUM = 0.                                                              
          DO j = 1, iY-1                                                        
           Efact = dabs(Eint3(TAU(iY)-TAU(j))-Eint3(TAU(iY)-TAU(j+1)))          
           fave = 0.5*(faux(j)+faux(j+1))                                       
           SUM = SUM + fave*Efact                                               
          END DO                                                                
          fp(iL,iY) = 0.5*SUM                                                   
c       and f(-) (in arg Tau>t)                                                 
          SUM = 0.                                                              
          DO j = iY, nY-1                                                       
           Efact = dabs(Eint3(TAU(iY)-TAU(j))-Eint3(TAU(iY)-TAU(j+1)))          
           fave = 0.5*(faux(j)+faux(j+1))                                       
           SUM = SUM + fave*Efact                                               
          END DO                                                                
          fm(iL,iY) = 0.5*SUM                                                   
        END DO                                                                  
                                                                                
        DO iY = 1, nY                                                           
           mat2(iL,iY) = fp(iL,iY) - fm(iL,iY)                                  
        END DO                                                                  
c     End of loop over iL                                                       
      END DO                                                                    

      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
      DOUBLE PRECISION FUNCTION Eint1(x)                                        
c     ======================================================================        
c     Needed for the slab geometry. It calculates the first exponential             
c     integral E1(x) by analytical f-la (13.13) from Abramovitz & Stegun(1994)      
c                                                             [MN,Dec'97]           
c     ======================================================================        
      IMPLICIT none                                                             
      INTEGER i                                                                 
      DOUBLE PRECISION AC(4),BC(4), CC(6), x, aux, poly, denom                  
      DATA AC/8.5733287401,18.0590169730,8.6347608925,0.2677737343/             
      DATA BC/9.5733223454,25.6329561486,21.0996530827,3.9584969228/            
      DATA CC/-0.57721566,0.99999193,-0.24991055,0.05519968,-0.00976004,        
     &        0.00107857/                                                       

c     For x=1E-15, E1~30 (used below to limit the value at x=0);for x>1, E1<1e-8   
c     Two approximations are used, for x>1 and x<1, respectively                   
      IF (x.GT.1.0) THEN                                                        
         poly = 0.0                                                             
         denom = 0.0                                                            
         aux = 1.0                                                              
         DO i = 1, 4                                                            
           poly = poly + AC(5-i)*aux                                            
           denom = denom + BC(5-i)*aux                                          
           aux = aux * x                                                        
         END DO                                                                 
           poly = poly + aux                                                    
           denom = denom + aux                                                  
         Eint1 = poly/denom/x*dexp(-x)                                          
      ELSE                                                                      
         IF (x.LE.1.E-15) x=1.0E-15                                             
         poly = 0.0                                                             
         aux = 1.0                                                              
         DO i = 1, 6                                                            
           poly = poly + CC(i)*aux                                              
           aux = aux * x                                                        
         END DO                                                                 
         Eint1 = poly - dlog(x)                                                 
      END IF                                                                    

      RETURN                                                                    
      END                                                                       
                                                                                

      DOUBLE PRECISION FUNCTION Eint2(x)                                        
c     ======================================================================        
c     Needed for the slab geometry. It calculates the second exponential            
c     integral E2(x) by the recurrence f-la. (see Abramovitz & Stegun,1994)         
c                                                             [MN,Dec'97]           
c     ======================================================================        
      IMPLICIT none                                                             
      DOUBLE PRECISION x, Eint1                                                 

       IF(x.LT.0.) x=dabs(x)                                                    
       IF (x.LT.1.0D-15) THEN                                                   
         Eint2 = 1.0                                                            
        ELSE                                                                    
         Eint2 = dexp(-x) - x*Eint1(x)                                          
       END IF                                                                   

      RETURN                                                                    
      END                                                                       

                                                                                
      DOUBLE PRECISION FUNCTION Eint3(x)                                        
c     ======================================================================        
c     Needed for the slab geometry. It calculates the third exponential             
c     integral E3(x) by the recurrence f-la. (see Abramovitz & Stegun,1994)         
c                                                            [MN,Dec'97]            
c     ======================================================================        
      IMPLICIT none                                                             
      DOUBLE PRECISION x, Eint1                                                 

       IF(x.LT.0.) x=dabs(x)                                                    
       IF (x.LT.1.0D-15) THEN                                                   
         Eint3 = 0.5                                                            
        ELSE                                                                    
         Eint3 = 0.5*((1.0-x)*dexp(-x)+x*x*Eint1(x))                            
       END IF                                                                   

      RETURN                                                                    
      END                                                                       
                                                                                

      DOUBLE PRECISION FUNCTION Eint4(x)                                        
c     ======================================================================        
c     Needed for the slab geometry. It calculates the fourth exponential            
c     integral E4(x) by the recurrence f-la. (see Abramovitz & Stegun,1994)         
c                                                             [MN,Jan'97]           
c     ======================================================================        
      IMPLICIT none                                                             
      DOUBLE PRECISION x, Eint3                                                 

       IF(x.LT.0.) x=dabs(x)                                                    
       IF (x.LT.1.0D-15) THEN                                                   
         Eint4 = 1.0d00/3.0d00                                                  
        ELSE                                                                    
         Eint4 = (dexp(-x)-x*Eint3(x))/3.0d00                                   
       END IF                                                                   

      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
      DOUBLE PRECISION FUNCTION Eint5(x)                                        
c     ======================================================================        
c     Needed for the slab geometry. It calculates the fifth exponential             
c     integral E5(x) by the recurrence f-la. (see Abramovitz & Stegun,1994)         
c                                                             [MN,Jan'97]           
c     ======================================================================        
      IMPLICIT none                                                             
      DOUBLE PRECISION x, Eint4                                                 

        IF(x.LT.0.) x=dabs(x)                                                   
       IF (x.LT.1.0d-15) THEN                                                   
         Eint5 = 0.25d00                                                        
        ELSE                                                                    
         Eint5 = 0.25d00*(dexp(-x)-x*Eint4(x))                                  
       END IF                                                                   

      RETURN                                                                    
      END                                                                       


      SUBROUTINE SLBacc(flux,accuracy,devmax,FbolOK,error)                      
c     =======================================================================       
c     This is SUB ChkFlux modified to be used for slab calculation.                 
c     Checks the bolometric flux conservation at any point of the slab grid.        
c     In case of nonconservation inserts a number of points at certain              
c     places.                                             [MN,99; ZI'96]            
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
      INTEGER iYins(npY), k, kins, iY, flag, error, istop, FbolOK               
      DOUBLE PRECISION Yins(npY), flux(npY), delTAUmax, devmax, ratio,          
     &        accuracy, devfac, ff, ffold                                       

      error = 0                                                                 
      kins = 0                                                                  
      devmax = 0.0                                                              
c     maximal delTAU is no more than 2 times the average value                  
      delTAUmax = 2.0*TAUtot(1)/nY                                              
c     find maximal relative deviation of fbol:                                  
      DO iY = 2, nY                                                             
        ratio = (dabs(flux(iY))-dabs(fmed))/dabs(fmed)                          
        IF (dabs(ratio).GT.devmax) devmax = dabs(ratio)                         
      END DO                                                                    
      ff = 0.0                                                                  
      istop = 0                                                                 
      devfac = 0.1                                                              
c     search for places to improve the grid                                     
      DO WHILE (istop.NE.1)                                                     
        DO iY = 2, nY                                                           
          ratio = (dabs(flux(iY))-dabs(fmed))/dabs(fmed)                        
          ffold = abs(flux(iY-1)/fmed-1.0)                                      
          ff = abs(flux(iY)/fmed-1.0)                                           
          flag = 0                                                              
c         if any of these criteria is satisfied insert a point:                 
          IF(ff.GT.accuracy) flag=1                                             
c         1) if error is increasing too fast                                    
          IF (abs(ff-ffold).GT.devfac*devmax) flag = 1                          
c         2) if delTAU is too large at the left edge:                           
          IF(TAUslb(1,iY).LT.5.) THEN                                           
           IF ((TAUslb(1,iY)-TAUslb(1,iY-1)).GT.                                
     &                             delTAUmax) flag = 1                          
          END IF                                                                
          IF(flag.EQ.1.AND.devmax.GE.accuracy) THEN                             
            kins = kins + 1                                                     
            Yins(kins) = Y(iY-1)+0.5*(Y(iY)-Y(iY-1))                            
            iYins(kins) = iY-1                                                  
          END IF                                                                
        END DO                                                                  
        IF (devmax.LT.accuracy.OR.devfac.LT.0.01) THEN                          
          istop = 1                                                             
        ELSE                                                                    
          IF (kins.GT.0) istop = 1                                              
        END IF                                                                  
        devfac = devfac / 2.0                                                   
      END DO                                                                    
                                                                                
      IF (kins.EQ.0) THEN                                                       
         FbolOK = 1                                                             
        ELSE                                                                    
c       Add all new points to Y(nY). This gives the new Y(nY+kins).             
c       However, check if npY is large enough to insert all points:             
        IF ((nY+kins).GT.npY) THEN                                              
         IF (iX.GE.1) THEN                                                      
         write(18,*)' ****************     WARNING   ******************'        
         write(18,*)'  The new Y grid can not accomodate more points!'          
         write(18,'(a,i5)')'   Specified accuracy would require',nY+kins        
         write(18,'(a,i5,a)')'   points, while npY =',npY,'.'                   
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

777   RETURN                                                                    
      END                                                                       
                                                                                

      SUBROUTINE SLBmisc(fbol,fmax,fmed,AveDev,RMS,nY)                          
c     =======================================================================       
c     Finds some additional quantities: fave-the average and fmed- the          
c     median value of fbol in the slab, average and RMS deviations of           
c     fbol from fmed.                                     [MN,Aug'98]           
c     =======================================================================       
      IMPLICIT NONE                                                             
      INTEGER iY, nY, imid, iup, idn                                            
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      DOUBLE PRECISION fbol(npY), fsort(npY), fave, fmed, AveDev,               
     &                 err, RMS, rnY, fmax, aux                                 
c     Find the median of fbol values distribution                               
      DO iY = 1, nY                                                             
        aux = DMAX1(dabs(fbol(iY)), accFbol)                                    
        fbol(iY) = aux                                                          
      END DO                                                                    
      DO iY = 1, nY                                                             
        fsort(iY) = fbol(iY)                                                    
      END DO                                                                    
      CALL sort(fsort,nY)                                                       
      IF (MOD(nY,2).NE.0) THEN                                                  
c       if nY odd:                                                              
        imid = (nY+1)/2                                                         
        fmed = fsort(imid)                                                      
      ELSE                                                                      
c       if nY even:                                                             
        iup = nY/2 + 1                                                          
        idn = nY/2                                                              
        fmed = 0.5*(fsort(iup) + fsort(idn))                                    
      END IF                                                                    
c     In sub PrOut: if(fmed.LE.accFbol) fmed =0                                 
      rnY = 1.0*nY                                                              
      fave = 0.                                                                 
      fmax = 0.                                                                 
c     find the max |fbol| (carried in slab.inc, needed in SLBacc)               
      DO iY = 1, nY                                                             
        fave = fave + fbol(iY)                                                  
        aux = dabs(fbol(iY))                                                    
        IF(aux.GT.fmax) fmax = aux                                              
      END DO                                                                    
      fave = fave / rnY                                                         
c     AveDev is the average relative deviation from the median                  
      AveDev = 0.                                                               
      RMS = 0.0                                                                 
      DO iY = 1, nY                                                             
        aux = dabs(fbol(iY))                                                    
c        err = dabs(fmed-aux) /fmax                                             
        err = dabs(fmed-aux) /fmed                                              
        AveDev = AveDev + err                                                   
        RMS = RMS + err*err                                                     
      END DO                                                                    
      RMS = sqrt(RMS/rnY/(rnY-1.))                                              
      AveDev = AveDev/rnY                                                       

      RETURN                                                                    
      END                                                                       

                                                                                
      SUBROUTINE SLBgray(model,nG,error)                                        
c     =======================================================================       
c     Solves the gray slab problem. For single lambda calculation               
c     comment the call to Spectral in main.for and the calls to FindErr         
c     in Analysis.for                                      [MN, May'98]         
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
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
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
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
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER iY, nG                                                            
      INTEGER error, nLst, model                                                
      character tt*58	                                                          
      DOUBLE PRECISION mat0(npL,npY,npY), Em(npL,npY), Uold(npL,npY),           
     &         fdsp(npL,npY),fdsm(npL,npY), fdep(npL,npY),fdem(npL,npY),        
     &         fp(npY), fm(npY), fd(npY)                                        

      nLst = nL                                                                 
      nL = 1                                                                    
c     generate some temporary arrays                                            
      DO iY = 1, nY                                                             
c       to prevent exp underflow:                                               
        IF(TAUslb(1,iY).GE.50.) THEN                                            
          Us(1,iY) = 0.                                                         
        ELSE                                                                    
          Us(1,iY) = dexp(-TAUslb(1,iY))                                        
        END IF                                                                  
        fs(1,iY) = Us(1,iY)                                                     
        Em(1,iY) = 0.                                                           
        omega(1,iY) = 1.                                                        
      END DO                                                                    
c     find radiative transfer matrices                                          
      IF (iX.GE.1) write(18,*)' Calculating matrix, pure scatt, iL=1'           
      CALL SLBMat(TAUslb,mat0)                                                  
c     solve for Utot                                                            
      CALL INVERT(nY,nL,mat0,Us,Em,omega,Utot,Uold,error)                       
c      find new Td                                                              
       DO iY = 1, nY                                                            
c       For a single lambda run:                                                
        Td(1,iY) = Tsub(1)*dsqrt(dsqrt(Utot(1,iY)/Utot(1,1)))                   
       END DO                                                                   
c      Find the diffuse scattered flux(fl=1 for scatt. and fl=0 is for emission)
       CALL SLBdifF(1,omega,TAUslb,Utot,nL,nY,fds,fdsp,fdsm)                    
c      Find the diffuse emitted flux                                            
       CALL SLBdifF(0,omega,TAUslb,Utot,nL,nY,fde,fdep,fdem)                    
       CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                                  
c      overall conservation of flux - min/max err                               
       DO iY = 1, nY                                                            
        fbol(iY) = fTot(1,iY)                                                   
        fsbol(iY) = fs(1,iY)                                                    
        fp(iY) = fdsp(1,iY)+fdep(1,iY)                                          
        fm(iY) = fdsm(1,iY)+fdem(1,iY)                                          
        fd(iY) = fp(iY)-fm(iY)                                                  
       END DO                                                                   
c     FindErr calculates the err acc. to min/max values of fbol                 
      CALL FindErr(fbol,maxFerr,nY)                                             
      CALL SLBmisc(fbol,fmax,fmed,AveDev,RMS,nY)                                

c     The slab Tau-profile at the fiducious lambda (needed in PrOut)            
      DO iY = 1, nY                                                             
        tr(iY) = TAUslb(1,iY)                                                   
      END DO                                                                    
      IF(iInn.EQ.1) THEN                                                        
       write(18,*) '----- Single Lambda Case (iL=1) -----'                      
       write(18,'(a10,1p,E11.3)')'    tauT =', TAUslb(1,nY)                     
        write(18,'(a11,1p,E11.3)')'     fmed =',fmed                            
        write(18,'(a11,1p,E11.3)')'  RMS err =',RMS                             
        write(18,'(a11,1p,E11.3)')'  MAX err =',maxFerr                         
                                                                                
       tt = '     tr       fTot(i)      fs(i)      fd(i)     fp(i)   '          
       write(18,'(a58,a30)') tt , '  fm(i)     Utot(i)      Td(i)'              
       DO iY = 1, nY                                                            
        write(18,'(1p,8E11.3)') tr(iY), fbol(iY), fsbol(iY), fd(iY),            
     &                          fp(iY), fm(iY), Utot(1,iY), Td(1,iY)            
       END DO                                                                   
      END IF                                                                    
c     this is for running many models with single lambda case                   
      SmC(5,model) = maxFerr                                                    
      nL = nLst                                                                 

      RETURN                                                                    
      END                                                                       


      SUBROUTINE SLBIniT(nG,alpha)                                              
c     =======================================================================       
c     This subroutine calculates the initial approximation for the temperature      
c     profile and alpha array. It is based on the analogous subroutine InitTemp     
c     for spherical shell written by ZI, Jul'96.                [MN, Dec'97]        
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
      DOUBLE PRECISION Sigext(npL), aux(npY), faux(npL), fstar(npY),            
     &       Qfstar(npY), fstarQ(npY), QP(npY), QF(npY), oldalpha(npY),         
     &       alpha(npG,npY), IntQF, xP, Planck, resaux, delta, FovrF1,          
     &       TAU(npY)                                                           

c     Lambda integral of fs=f_e*exp(-TAU) -> fstar                              
      CALL Bolom(fs,fstar)                                                      
c     loop over grains                                                          
      DO iG = 1, nG                                                             
c       first approximation for temperature                                     
        DO iY = 1, nY                                                           
          Td(iG,iY) = Tsub(iG)                                                  
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
            faux(iL) = Sigext(iL) *fs(iL,iY)/lambda (iL)                        
          END DO                                                                
          CALL Simpson(npL,1,nL,lambda,faux,resaux)                             
          Qfstar(iY) = resaux                                                   
        END DO                                                                  
c      The slab geom. needs initial approximation for F/F1; Psi is alpha(nG,1)  
c      (The solution is not sensitive to this approximation, MN)                
        IF(TAUmax.LT.500.) THEN                                                 
           FovrF1 = 1.0                                                         
         ELSE                                                                   
           FovrF1 = 1. - 0.25*alpha(1,1)                                        
        END IF                                                                  
c       iterate until Td and Psi converge (i.e. until alpha does)               
        istop = 0                                                               
        DO WHILE (istop.NE.1)                                                   
          istop = 1                                                             
          DO iY = 1, nY                                                         
c           generate auxiliary function for lambda integration:                 
            DO iL = 1, nL                                                       
              xP = 14400.0 / Td(iG,iY) / lambda(iL)                             
              faux(iL) = Sigext(iL) * Planck(xP) / lambda (iL)                  
            END DO                                                              
            CALL Simpson(npL,1,nL,lambda,faux,resaux)                           
c           Planck average of Sigext                                            
            QP(iY) = resaux                                                     
c           calculate QF for slab (IV in the notes, the equivalent of eq.B5)    
            QF(iY) = Qfstar(iY)/FovrF1+QP(iY)*(1.0-fstar(iY)/FovrF1)            
c           Find the second term in III (in the notes)                          
            DO iL = 1, nL                                                       
              faux(iL) = fs(iL,iY)*(1-SigmaA(iG,iL)/QP(iY))/lambda(iL)          
            END DO                                                              
            CALL Simpson(npL,1,nL,lambda,faux,resaux)                           
            fstarQ(iY) = resaux                                                 
c           store current alpha                                                 
            oldalpha(iY) = alpha(iG,iY)                                         
          END DO                                                                
c         for each y calculate new alpha from III (the equivalent of eq.B7)     
          DO iY = 1, nY                                                         
            TAU(iY) = TAUslb(iLfid,iY)                                          
          END DO                                                                
          DO iY = 1, nY                                                         
            DO iW = iY, nY                                                      
              aux(iW) = FovrF1*QF(iY)                                           
            END DO                                                              
            CALL Simpson(npY,iY,nY,TAU,aux,IntQF)                               
            alpha(iG,iY) = 3.*IntQF - fstarQ(iY)                                
            alpha(iG,iY) = oldalpha(iY)                                         
c           calculate temperature                                               
            Td(iG,iY) = Tsub(iG) * (alpha(iG,iY)/alpha(iG,1))**0.25             
            delta = DABS((alpha(iG,iY)-oldalpha(iY))/alpha(iG,iY))              
            IF (delta.GT.accConv) istop = 0                                     
          END DO                                                                
c       end of iterations                                                       
        END DO                                                                  
c     end of loop over grains (iG)                                              
      END DO                                                                    

      RETURN                                                                    
      END                                                                       
                                                                                

      SUBROUTINE SLBmat(TAUslb,mat0)                                            
c     ======================================================================        
c     This subroutine evaluates the radiative transfer matrix for slab geometry.    
c     TAUslb is the array of optical depths along the line of sight. [MN, Dec'97]   
c     ======================================================================        
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
      INTEGER iL, iY, i, j, k                                                   
      DOUBLE PRECISION TAUslb(npL,npY), TAUr(npY), mat0(npL,npY,npY),           
     &         L1(npY,npY),L2(npY,npY),L3(npY,npY),L4(npY,npY),                 
     &         alpha(npY,npY), beta(npY,npY), gamma(npY,npY),                   
     &         delta(npY,npY), SUM                                              

      DO iL = 1, nL                                                             
       DO j = 1, nY                                                             
        DO k = 1, nY                                                            
          mat0(iL,j,k) = 0.0                                                    
        END DO                                                                  
       END DO                                                                   
      END DO                                                                    
c     -- evaluate matrix elements --                                            
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c      'Myspline' needs a 1D array                                              
        DO iY = 1, nY                                                           
          TAUr(iY) = TAUslb(iL,iY)                                              
        END DO                                                                  
c       Get the spline coefficients of the source function                      
        CALL MYSPLINE(TAUr,nY,alpha,beta,gamma,delta)                           
        CALL SLBAnalyt(TAUr,L1,L2,L3,L4,nY,npY)                                 
c       Matrix element for U:                                                   
        DO k = 1, nY                                                            
          DO j = 1, nY                                                          
            SUM = 0.                                                            
            DO i = 1, nY                                                        
             SUM = SUM + L1(i,k)*alpha(i,j)+L2(i,k)*beta(i,j)                   
     &                 + L3(i,k)*gamma(i,j)+L4(i,k)*delta(i,j)                  
            END DO                                                              
            mat0(iL,k,j) = SUM                                                  
          END DO                                                                
        END DO                                                                  
c     end of loop over wavelengths                                              
      END DO                                                                    
c     save Y to Yok; needed for analysis in cases when requirement for          
c     finer grids cannot be satisfied and previous solution is used for         
c     output                                                                    
      nYok = nY                                                                 
      DO iY = 1, nY                                                             
        Yok(iY) = Y(iY)                                                         
      END DO                                                                    

      RETURN                                                                    
      END                                                                       

                                                                                
      SUBROUTINE SLBRadT(nG,error)                                              
c     ======================================================================        
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
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
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
      INTEGER iY, nG, itnum, itlim                                              
      INTEGER error, Conv, iter, Fconv, Uconv, BolConv, m14                     
      DOUBLE PRECISION mat0(npL,npY,npY), Em(npL,npY), alpha(npG,npY),          
     &    fdsp(npL,npY),fdsm(npL,npY),fdep(npL,npY),fdem(npL,npY),dmaxU,        
     &    fbolold(npY),Uold(npL,npY), dmaxF, olderr, oldmaxU                    

      IF (iInn.eq.1) THEN                                                       
        write(38,'(a8,i5)') '    nY= ',nY                                       
        write(38,*) '    iter   maxFerr     dmaxU       dmaxF'                  
      END IF                                                                    
c     The slab Tau-profile at the fiducious lambda (needed in PrOut)            
      DO iY = 1, nY                                                             
        tr(iY) = TAUslb(iLfid,iY)/TAUslb(iLfid,nY)                              
      END DO                                                                    
c     generate stellar moments                                                  
      CALL SLBStar(error)                                                       
c     generate albedo through the envelope                                      
      CALL getOmega(nG)                                                         
c     finish when file with the stellar spectrum is not available               
      IF (error.EQ.3) goto 999                                                  
c     generate the first approximations for Td and alpha                        
      CALL SLBiniT(nG,alpha)                                                    
c     find radiative transfer matrices                                          
      IF (iX.GE.1) write(18,*)' Calculating weight matrices'                    
      CALL SLBmat(TAUslb,mat0)                                                  
      IF (iX.GE.1) THEN                                                         
        write(18,*)' Weight matrices OK, calculating Tdust'                     
      END IF                                                                    
      itlim = 9000                                                              
      itnum = 4                                                                 
      Conv = 0                                                                  
      iter = 0                                                                  
      oldmaxU = 0.                                                              
      olderr = 0.                                                               
      m14 = 0                                                                   
c     === Iterations over dust temperature =========                            
      DO WHILE (Conv.EQ.0.AND.iter.LE.itlim)                                    
        iter = iter + 1                                                         
c       find emission term                                                      
        CALL Emission(0,nG,Td,alpha,abund,Us,Em)                                
c       solve for Utot                                                          
        CALL Invert(nY,nL,mat0,Us,Em,omega,Utot,Uold,error)                     
c       find new Td and alpha, and check convergence                            
        CALL FindTemp(0,Utot,nG,Td,alpha)                                       
c       --------------------------------------                                  
c       every itnum-th iteration check convergence:                             
c       first find 'old' flux (i.e. in the previous iteration)                  
        IF (MOD(iter+1,itnum).EQ.0) THEN                                        
c         Find the diffuse scattered flux(fl=1 for scatt. and fl=0 is for emissi
          CALL SLBdifF(1,omega,TAUslb,Utot,nL,nY,fds,fdsp,fdsm)                 
c         Find the diffuse emitted flux                                         
          CALL SLBdifF(0,omega,TAUslb,Em,nL,nY,fde,fdep,fdem)                   
          CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                               
c         find bolometric flux                                                  
          CALL Bolom(ftot,fbolold)                                              
        END IF                                                                  
        IF (MOD(iter,itnum).EQ.0) THEN                                          
c         Find the diffuse scattered flux                                       
          CALL SLBdifF(1,omega,TAUslb,Utot,nL,nY,fds,fdsp,fdsm)                 
c         Find the diffuse emitted flux                                         
          CALL SLBdifF(0,omega,TAUslb,Em,nL,nY,fde,fdep,fdem)                   
c         Add them to the stellar flux to find total flux                       
          CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                               
c         find bolometric flux                                                  
          CALL Bolom(ftot,fbol)                                                 
c         check convergence of bolometric flux                                  
          CALL Converg1(nY,accFbol,dynrange,fbolold,fbol,Fconv,dmaxF)           
c         check convergence of energy density                                   
          CALL Converg2(nY,nL,accConv,dynrange,Uold,Utot,Uconv,dmaxU)           
c         FindErr calculates the err acc. to min/max values                     
          CALL FindErr(fbol,maxFerr,nY)                                         
c         flag m14 is to prevent unnecessary looping if convergence             
c         is extremely slow and err does not improve.                           
          IF(iter.GE.1000.AND.                                                  
     &           dabs(dmaxU-oldmaxU).LE.1.0d-5.AND.                             
     &             dabs(maxFerr-olderr).LE.1.0d-4) THEN                         
            m14 = 1                                                             
            Conv = 1                                                            
          END IF                                                                
          olderr = maxFerr                                                      
          oldmaxU = dmaxU                                                       
c         ------  printout of errors and convergence with iter.(inner flag): ----       
          IF(iInn.EQ.1) THEN                                                    
            write(38,'(i7,1p,3e12.4)') iter, maxFerr, dmaxU, dmaxF              
          END IF                                                                
c         ------------------------------------------------------------------------      
          IF (maxFerr.LE.accuracy) THEN                                         
            BolConv = 1                                                         
          ELSE                                                                  
            BolConv = 0                                                         
          END IF                                                                
c         total criterion for convergence: Utot must converge, and ftot         
c         must either converge or have the required accuracy                    
          IF (Uconv*(Fconv+BolConv).GT.0) Conv = 1                              
        END IF                                                                  
      END DO                                                                    
c     === The End of Iterations over Td ===                                     
      IF (iX.GE.1) THEN                                                         
        IF (iter.le.itlim) THEN                                                 
          write(18,*)' Convergence achieved, number of'                         
          write(18,'(a34,i4)')                                                  
     &      ' iterations over energy density: ',iter                            
        ELSE                                                                    
          write(18,'(a38,i5,a6)')                                               
     &    '  No convergence on energy density in',iter,' iter.'                 
          iWARNING = iWARNING + 1                                               
        END IF                                                                  
        write(18,'(a25,1p,E9.2)')                                               
     &    '  Max error in bol.flux:',maxFerr                                    
      END IF                                                                    
c     calculate the emission term for the converged Td                          
      CALL Emission(0,nG,Td,alpha,abund,Us,Em)                                  
c     calculate flux                                                            
      CALL SLBdifF(1,omega,TAUslb,Utot,nL,nY,fds,fdsp,fdsm)                     
      CALL SLBdifF(0,omega,TAUslb,Em,nL,nY,fde,fdep,fdem)                       
      CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                                   
      CALL Add2(fdsp,fdep,fpbol,nY)                                             
      CALL Add2(fdsm,fdem,fmbol,nY)                                             
      CALL Bolom(ftot,fbol)                                                     
      CALL Bolom(fs,fsbol)                                                      
c     Find the err acc. to min/max values of fbol:                              
      CALL FindErr(fbol,maxFerr,nY)                                             
c     find average value of fbol and some misceleneous quantities               
      CALL SLBmisc(fbol,fmax,fmed,AveDev,RMS,nY)                                
c     calculate additional output quantities                                    
      CALL Multiply(1,npY,nY,npL,nL,mat0,Utot,omega,0,Us,Uds)                   
      CALL Multiply(0,npY,nY,npL,nL,mat0,Em,omega,0,Us,Ude)                     
      CALL Add(npY,nY,npL,nL,Us,Uds,Ude,Uchck)                                  
      CALL Bolom(Utot,Ubol)                                                     
      CALL Bolom(Uchck,UbolChck)                                                
c     fmed is printed in .out file as f1 = fmed = F/Fe1                         
c     ============ if the inner flag iInn=1:  =========                             
      IF(iX.GE.1 .AND. iInn.EQ.1) THEN                                          
        write(18,'(a11,1p,E11.3)')'   TAUfid =',TAUfid                          
        write(18,'(a11,1p,E11.3)')'     fmed =',fmed                            
        write(18,'(a11,1p,E11.3)')'     fmax =',fmax                            
        write(18,'(a11,1p,E11.3)')'  MAX err =',maxFerr                         
        write(18,*)                                                             
     &'     tr      fbol       fsbol      fdbol       Ubol     fpbol'           
        DO iY = 1, nY                                                           
         write(18,'(1p,6E11.3)') tr(iY), fbol(iY), fsbol(iY),                   
     &               (fpbol(iY)-fmbol(iY)),Ubol(iY),fpbol(iY)                   
        END DO                                                                  
      END IF                                                                    
c     ===================================================                           

999   RETURN                                                                    
      END                                                                       


      SUBROUTINE SLBsolve(model,nG,error)                                       
c     =======================================================================       
c     This subroutine solves the continuum radiative transfer problem in            
c     planar geometry.                                      [MN, Dec.'97]           
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
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER model, error, nG, iterFbol, FbolOK, grid                          
      DOUBLE PRECISION devmax, oldFerr, mu, h, k                                

c     counter for iterations over bolometric flux conservation                  
      iterFbol = 0                                                              
      FbolOK = 0                                                                
      oldFerr = 100.                                                            
c     make the grid for the smaller mu                                          
      IF (dabs(mu1).LE.dabs(mu2)) THEN                                          
        mu = dabs(mu1)                                                          
      ELSE                                                                      
        mu = dabs(mu2)                                                          
      END IF                                                                    
      IF((mu.LT.0.1).AND.(TAUmax.GE.20.)) THEN                                  
        grid = 1                                                                
c      oblique illumination grid (grid=1) with initial parameters               
       IF(npY.LT.50) THEN                                                       
        iWARNING = iWARNING + 1                                                 
        error = 4                                                               
        RETURN                                                                  
       ELSE                                                                     
        nY = 50                                                                 
       END IF                                                                   
        h = 0.5                                                                 
        k = 5.0                                                                 
        CALL SLBtrace(TAUmax,TAUTot,TAUslb,Y,mu,h,k,nL,nY)                      
      ELSE                                                                      
        grid = 0                                                                
c       generate the increment for the initial TAU-grid (grid=0)                
        CALL SLBY(TAUmax,Y,nY)                                                  
      END IF                                                                    
c     ------------- Loop over bol.flux conservation -------------                   
      DO WHILE (FbolOK.EQ.0)                                                    
       iterFbol = iterFbol + 1                                                  
        IF((iterFbol.GE.4).AND.(grid.EQ.1).AND.(iX.GE.1)) THEN                  
          write(18,*)'  Can not improve the oblique'                            
          write(18,*)'  illumination grid anymore  '                            
          iWARNING = iWARNING + 1                                               
          FbolOK = 2                                                            
          RETURN                                                                
        END IF                                                                  
        IF (iX.GE.1) THEN                                                       
          write(18,*)'  ',iterFbol,' iteration over Fbol'                       
        END IF                                                                  
        IF (iVerb.EQ.2) write(*,'(a14,i3,a20)')                                 
     &     ' In SLBsolve: ',iterFbol,' iteration over Fbol'                     
c       if the illumination angle is too large choose the proper grid           
        IF(grid.NE.1) THEN                                                      
          CALL SLBtau(TAUTot,TAUslb,Y,nL,nY)                                    
        END IF                                                                  
        IF (iX.GE.1) THEN                                                       
          write(18,'(a28,i3,a12)') '  Calculation for slab with ',              
     &                          nY,' grid points'                               
        END IF                                                                  
        IF (nY.GT.npY) THEN                                                     
         IF (iVerb.EQ.2) write(*,*)' npY needs to be increased!'                
         IF (iX.GE.1) THEN                                                      
          write(18,*) ' ********** MESSAGE from SOLVE *********'                
          write(18,*) ' npY has to be increased. Please, choose'                
          write(18,*) ' parameters for slab in file userpar.inc'                
          write(18,*) ' ***************************************'                
          iWARNING = iWARNING + 1                                               
          RETURN                                                                
         END IF                                                                 
        END IF                                                                  
c       solve the radiative transfer problem                                    
        CALL SLBRadT(nG,error)                                                  
         IF (iVerb.EQ.2) write(*,*)'Done with SLBRadT'                          
c       for tests of single lambda case:                                        
c       CALL SLBgray(model,nG,error)                                            
        IF(maxFerr.LE.accuracy) THEN                                            
          FbolOK = 1                                                            
        ELSE                                                                    
         IF(grid.EQ.1) THEN                                                     
           IF(iterFbol.EQ.1) THEN                                               
c            next try with 104pt grid:                                          
             h = h/1.5                                                          
             k = k*1.5                                                          
             nY = 104                                                           
             CALL SLBtrace(TAUmax,TAUTot,TAUslb,Y,mu,h,k,nL,nY)                 
           END IF                                                               
           IF(iterFbol.EQ.2) THEN                                               
c            last attempt with 160pt grid:                                      
             h = 0.25                                                           
             k = 7.5                                                            
             nY = 160                                                           
             CALL SLBtrace(TAUmax,TAUTot,TAUslb,Y,mu,h,k,nL,nY)                 
           END IF                                                               
         ELSE                                                                   
           CALL SLBacc(fbol,accuracy,devmax,FbolOK,error)                       
           IF(iterFbol.GT.1.AND.maxFerr.GT.oldFerr) THEN                        
            IF(iX.GE.1) THEN                                                    
             write(18,*)' ================ WARNING ================= '          
             write(18,*)' Error is increasing in spite of grid'                 
             write(18,*)' improvement. Stopping iterations over Fbol.'          
            END IF                                                              
            error = 0                                                           
            FbolOK = 2                                                          
            iWARNING = iWARNING + 1                                             
           END IF                                                               
           oldFerr = maxFerr                                                    
         END IF                                                                 
c         if the grid size limit is reached error=2                             
          IF (error.EQ.2.AND.iterFbol.EQ.1) THEN                                
c          if this is the first calculation end this model                      
           IF(iX.GE.1) THEN                                                     
            write(18,*)' =========== IMPORTANT WARNING =========== '            
            write(18,*)' The limit for grid size is already reached'            
            write(18,*)' and flux conservation can not be improved'             
            write(18,'(a,1p,e9.2)')'  Max deviation of Fbol ',maxFerr           
            write(18,*)' Treat all results with caution!'                       
           END IF                                                               
           error = 0                                                            
           FbolOK = 2                                                           
           iWARNING = iWARNING + 1                                              
          END IF                                                                
c         if this is a higher iteration use previous solution                   
          IF (error.EQ.2) THEN                                                  
            IF (iX.GE.1.AND.iterFbol.GT.1) THEN                                 
              write(18,*)' ======= IMPORTANT WARNING ======== '                 
              write(18,*)' In trying to conserve Fbol reached'                  
              write(18,*)' the limit for grid sizes.  '                         
              write(18,'(a,1p,e9.2)')'  Max deviation of Fbol:',maxFerr         
              write(18,*)' Treat all results with caution!'                     
            END IF                                                              
            error = 0                                                           
            FbolOK = 2                                                          
            iWARNING = iWARNING + 1                                             
          END IF                                                                
c         if Fbol not conserved try again with a finer grid                     
          IF (FbolOK.EQ.0 .AND. iX.GE.1) THEN                                   
            write(18,*)'  ******** MESSAGE from SOLVE ********'                 
            write(18,'(a,1p,e9.2)')                                             
     &                '  Max deviation of Fbol:', maxFerr                       
            write(18,*)'  Trying again with finer grids'                        
          END IF                                                                
c        if could not conserve Fbol in 10 trials give it up                     
         IF (FbolOK.EQ.0 .AND. iterFbol.GE.10) THEN                             
          IF (iX.GE.1) THEN                                                     
          write(18,*)' **********  WARNING from SOLVE  **********'              
          write(18,*)' Could not obtain required accuracy in 10 trials.'        
          write(18,'(a,1p,e9.2)')'  Max deviation of Fbol:',maxFerr             
          write(18,*)' !!!!  Treat all results with caution  !!!!'              
          write(18,*)' ******************************************'              
          END IF                                                                
          iWARNING = iWARNING + 1                                               
          FbolOK = 2                                                            
         END IF                                                                 
        END IF                                                                  
c     end of loop over flux conservation                                        
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE SLBStar(error)                                                 
c     =======================================================================       
c     This subroutine generates the stellar moments in case of slab geometry.       
c                                                              [MN, Feb.'99]        
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
      INTEGER iY,iL,ios1,iLs,nLs, k,kstop,i,is, error,isn, Nlambdam             
c     Nlambdam is the max number entries for a user supplied stellar spectrum   
      PARAMETER (Nlambdam = 10000)	                       	                     
      DOUBLE PRECISION lambdaS(Nlambdam), Llamstar(Nlambdam), Lstar,            
     &       Stellar(Nlambdam), llS(Nlambdam), lS(Nlambdam), fsrc(npL),         
     &       fpl(npL), Bb(2), arg, EMfunc, a, b, x, Planck, BP,                 
     &       dyn2, UsL, UsR, Eint2, Eint3                                       

      IF(ksi.GT.0) THEN                                                         
        isn = 2                                                                 
      ELSE                                                                      
        isn = 1                                                                 
      END IF                                                                    
      dyn2 = dynrange*dynrange                                                  
      DO is = 1, isn                                                            
c      if startyp.GE.4 stellar spectrum is read from file 'nameStar'            
       IF (startyp(is).GE.4.AND.startyp(is).LE.6) THEN                          
        open(3,ERR=998,file=nameStar(is),STATUS='OLD')                          
        rewind(3)                                                               
        read(3,'(a235)',ERR=998) line                                           
        read(3,'(a235)',ERR=998) line                                           
        read(3,'(a235)',ERR=998) line                                           
        ios1 = 0                                                                
        iLs = 0                                                                 
        DO WHILE (ios1.ge.0)                                                    
          read(3,*,END=900,ERR=998,iostat=ios1) a,b                             
          IF (ios1.ge.0) THEN                                                   
            iLs = iLs + 1                                                       
            lambdaS(iLs) = a                                                    
            IF (a.LE.0.0) goto 998                                              
c           it is assumed that Llamstar is L_Lambda, but...                     
c           if startyp.EQ.3 then file gives lambda*L_lambda                     
            IF (startyp(is).EQ.4) Llamstar(iLs) = b / a                         
c           if startyp.EQ.4 then file gives L_lambda                            
            IF (startyp(is).EQ.5) Llamstar(iLs) = b                             
c           if startyp.EQ.5 then file gives Lnu=lambda**2*L_lambda              
            IF (startyp(is).EQ.6) Llamstar(iLs) = b / a / a                     
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
        IF (startyp(is).EQ.3) THEN                                              
          fsrc(is) = 1.0                                                        
          IF (Nlamtr(is).GT.1) THEN                                             
            DO i = 2, Nlamtr(is)                                                
              fsrc(i) = fsrc(i-1)*                                              
     &                        (lamtr(is,i-1)/lamtr(is,i))**klam(is,i-1)         
            END DO                                                              
          END IF                                                                
          DO iL = 1, nL                                                         
            IF((lambda(iL)-lamtr(is,1))*(lambda(iL)-                            
     &                             lamtr(is,Nlamtr(is)+1)).LE.0.0) THEN         
              kstop = 0                                                         
              k = 0                                                             
              DO WHILE (kstop.EQ.0)                                             
                k = k + 1                                                       
                IF (lambda(iL).GE.lamtr(is,k)) THEN                             
                  kstop = 1                                                     
                  fpl(iL) = fsrc(k)*(lamtr(is,k)/lambda(iL))**klam(is,k)        
                END IF                                                          
              END DO                                                            
            ELSE                                                                
              fpl(iL) = 0.0                                                     
            END IF                                                              
          END DO                                                                
        END IF                                                                  
       END IF                                                                   
       DO iY = 1, nY                                                            
c       loop over wavelengths                                                   
        DO iL = 1, nL                                                           
         IF (startyp(is).EQ.1) THEN                                             
            Bb(is) = 0.0                                                        
            DO k = 1, nBB(is)                                                   
              x = 14400.0 / lambda(iL) / Tbb(is,k)                              
              Bb(is) = Bb(is) + rellum(is,k)*Planck(x)                          
            END DO                                                              
         ELSE IF (startyp(is).EQ.2) THEN                                        
            Bb(is) = EMfunc(lambda(iL),Tbb(1,1),xSiO)                           
         ELSE IF (startyp(is).EQ.3) THEN                                        
            Bb(is) = fpl(iL)                                                    
         ELSE IF (lambda(iL).GT.lambdaS(nLs)) THEN                              
c           for lambda longer than the longest entry in nameStar                
c           assume Rayleigh-Jeans tail                                          
            Bb(is) = Stellar(nLs)*(lambdaS(nLs)/lambda(iL))**3.                 
         ELSE IF (lambda(iL).LT.lambdaS(1)) THEN                                
c           if shorter than the shortest assume 0                               
            Bb(is) = 0.0                                                        
         ELSE                                                                   
c          from a file: if within limits interpolate                            
           CALL LinInter(Nlambdam,nLs,lambdaS,Stellar,lambda(iL),iLs,BP)        
           Bb(is) = BP                                                          
         END IF                                                                 
c        ----------  done with stellar spectrum ---------------                      
c        stellar part of flux and en.density                                    
         IF (is.EQ.1) THEN                                                      
           IF((mu1+1.).LE.1.0e-5) THEN                                          
             fsL(iL,iY) = Bb(1)*Eint3(TAUslb(iL,iY))                            
           ELSE                                                                 
             x = TAUslb(iL,iY) / mu1                                            
             IF(x.GE.50.) THEN                                                  
               fsL(iL,iY) = 0.                                                  
             ELSE                                                               
               fsL(iL,iY) = Bb(1)*exp(-x)                                       
             END IF                                                             
           END IF                                                               
         END IF                                                                 
         IF (is.EQ.2) THEN                                                      
           IF((mu2+1.).LE.1.0e-5) THEN                                          
            arg = TAUslb(iL,nY)-TAUslb(iL,iY)                                   
            fsR(iL,iY) = Bb(2)*Eint3(arg)                                       
           ELSE                                                                 
             x = (TAUslb(iL,nY)-TAUslb(iL,iY)) / mu2                            
             IF(x.GE.50.) THEN                                                  
                 fsR(iL,iY) = 0.                                                
             ELSE                                                               
               fsR(iL,iY) = Bb(2)*exp(-x)                                       
             END IF                                                             
           END IF                                                               
         END IF                                                                 
        END DO                                                                  
       END DO                                                                   
c     end do over 'is' - the counter for sources                                
      END DO                                                                    
                                                                                
      DO iY = 1, nY                                                             
c       loop over wavelengths                                                   
        DO iL = 1, nL                                                           
          fs(iL,iY) = fsL(iL,iY) - ksi * fsR(iL,iY)                             
c         find en.density if diffuse illumination on the left (mu1=-1)          
          IF((mu1+1.).LE.1.0e-5) THEN                                           
c           to prevent exp. underflow                                           
            IF (TAUslb(iL,iY).GE.300.) THEN                                     
             UsL = 0.                                                           
            ELSE                                                                
             UsL = fsL(iL,iY)/Eint3(TAUslb(iL,iY))*Eint2(TAUslb(iL,iY))         
            END IF                                                              
          ELSE                                                                  
            UsL = fsL(iL,iY)/mu1                                                
          END IF                                                                
c         if diffuse illumination from the right                                
          IF((mu2+1.).LE.1.0e-5) THEN                                           
            arg = TAUslb(iL,nY)-TAUslb(iL,iY)                                   
c           to prevent exp. underflow                                           
            IF (arg.GE.300.) THEN                                               
             UsR = 0.                                                           
            ELSE                                                                
             UsR = fsR(iL,iY)/Eint3(arg)*Eint2(arg)                             
            END IF                                                              
          ELSE                                                                  
            UsR = fsR(iL,iY)/mu2                                                
          END IF                                                                
           Us(iL,iY) = UsL + ksi * UsR                                          
c         here only Us needs limit from below; fs can be negative though!       
          IF (Us(iL,iY).LT.dyn2) Us(iL,iY) = 0.0                                
        END DO                                                                  
      END DO                                                                    
                                                                                
c     normalize stellar quantities with the stellar bolometric flux             
      CALL Bolom(fsL,fsbol)                                                     
c     fsbol(1) is always non-zero                                               
      DO iY = 1, nY                                                             
        DO iL = 1, nL                                                           
           Us(iL,iY) = Us(iL,iY) / fsbol(1)                                     
           fs(iL,iY) = fs(iL,iY) / fsbol(1)                                     
           fsL(iL,iY) = fsL(iL,iY) / fsbol(1)                                   
           fsR(iL,iY) = fsR(iL,iY) / fsbol(1)                                   
        END DO                                                                  
      END DO                                                                    
      error = 0                                                                 
      goto 999                                                                  
998   write(12,*)' *** FATAL ERROR IN DUSTY! *************************'         
      write(12,*)' File with the spectral shape of external radiation:'         
      write(12,'(a2,a70)')'  ',nameStar(is)                                     
      write(12,*)' is missing or not properly formatted?!'                      
      write(12,*)' ***************************************************'         
      error = 3                                                                 

999   RETURN                                                                    
      END                                                                       

                                                                                
      SUBROUTINE SLBTau(TauTot,TAUslb,delT,nL,nY)                               
c     ======================================================================        
c     It generates TAUslb(npL,npY) grid with given spacing delT, calculated         
c     in sub SLBy. In the current version it is exp near the slab faces and         
c     equidistant in the middle.			             [MN,Sep'98]                         
c     ======================================================================        
      IMPLICIT none                                                             
      INTEGER nL, nY, iY, iL                                                    
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      DOUBLE PRECISION TAUslb(npL,npY),TauTot(npL),TAU(npY),delT(npY)           

       TAU(1) = 0.                                                              
       DO iY = 2, nY                                                            
         TAU(iY) = TAU(iY-1) + delT(iY)                                         
       END DO                                                                   
      DO iL = 1, nL                                                             
        DO iY = 1, nY                                                           
          TAUslb(iL,iY) = TauTot(iL)*TAU(iY)/TAU(nY)                            
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       

      SUBROUTINE SLBtrace(TAUmax,TAUTot,TAUslb,Y,mu,h,k,nL,nY)                  
c     ======================================================================        
c         This subroutine generates the TAU-grid for slab in the case of            
c         grazing angle of incidence.                            [MN,Jan'99]        
c     ======================================================================        
      IMPLICIT NONE                                                             
      INTEGER nL, nY, iY, iL, Nthick, Nskin, iadj                               
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      DOUBLE PRECISION TAUslb(npL,npY), TAUTot(npL), TAU(npY), Y(npY),          
     &              fsort(npY), TAUmax, mu, h, k, step                          

c     Number of pts. in the skin layer with depth=min{mu1,mu2}:                 
      Nskin = k/h                                                               
c     Counter for the end of the layer adjacent to the skin layer               
      iadj = 2*Nskin                                                            
c     Number of pts. in the opt.thick inner part:                               
      Nthick = nY - 4*Nskin                                                     
      iY = 1                                                                    
      TAU(1) = 0.0                                                              
      TAU(nY) = 1.0                                                             
c     resolve the surface boundary layer, where tau = 0 .. k*mu                 
      DO WHILE(iY.LT.nY/2)                                                      
      iY = iY + 1                                                               
        IF((TAU(iY).LE.k*mu).AND.(iY.LE.Nskin)) THEN                            
          TAU(iY) = TAU(1) + (iY-1)*h*mu/TAUmax                                 
          TAU(nY-iY+1) = TAU(nY) - TAU(iY)                                      
        ELSE                                                                    
c         resolve the adjacent layer                                            
          IF(iY.LE.iadj) THEN                                                   
            TAU(iY) = TAU(iY-1) + h/TAUmax                                      
            TAU(nY-iY+1) = TAU(nY) - TAU(iY)                                    
          ELSE                                                                  
c           for the opt.thick inner part:                                       
            step = (TAU(nY-iadj+1)-TAU(iadj)) / (Nthick+1)                      
            TAU(iY) = TAU(iY-1) + step                                          
            TAU(nY-iY+1) = TAU(nY) - TAU(iY)                                    
          END IF                                                                
        END IF                                                                  
      END DO                                                                    
c     sorting (to avoid overlapping or mismatched pieces)                       
      DO iY = 1, nY                                                             
        fsort(iY) = TAU(iY)                                                     
      END DO                                                                    
      CALL sort(fsort,nY)                                                       
      DO iY = 1, nY                                                             
        TAU(iY) = fsort(iY)                                                     
        Y(iY) = TAU(iY)                                                         
      END DO                                                                    
      DO iL = 1, nL                                                             
        DO iY = 1, nY                                                           
          TAUslb(iL,iY) = TAUTot(iL)*TAU(iY)                                    
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE SLBy(TAUmax,dTAU,nY)                                           
c     =======================================================================       
c     This subroutine generates the increment for the initial TAU-grid.         
c     It is equidistant in TAU in the middle and equidistant in log(TAU)        
c     near the two faces. It has 15 pts for small (tauV.LE.1) and 30pts         
c     for large total opt.depth (tauV.LE.1000). In some cases of tauV.GE.500.   
c     80pts grid is better to start with. If this is the case, comment          
c     the larger value of lim and uncomment lim=500, which will switch to       
c     80pt grid. Be sure you have selected the line for slab calculation        
c     in 'userpar.inc'.                                      [MN,Sep'98]        
c     =======================================================================       
      IMPLICIT NONE                                                             
      INTEGER nY, iY                                                            
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      DOUBLE PRECISION dTAU(npY), TAUmax                                        

      IF(TAUmax.LE.10.) THEN                                                    
          nY = 15                                                               
          DO iY = 1, nY                                                         
           IF(iY.LE.5) THEN                                                     
             dTAU(iY) = 0.5*(exp(0.5*iY)-exp(0.5))                              
           ELSE IF(iY.GT.5.and.iY.LE.12) THEN                                   
             dTAU(iY) = dTAU(iY-1)                                              
           ELSE                                                                 
             dTAU(iY) = dTAU(2+nY-iY)                                           
           END IF                                                               
          END DO                                                                
      ELSE                                                                      
       IF(TAUmax.LE.10000.) THEN                                                
         nY = 30                                                                
         DO iY = 1, nY                                                          
           IF(iY.LE.10) THEN                                                    
             dTAU(iY) = 0.5*(exp(0.5*iY)-exp(0.5))                              
            ELSE IF(iY.GT.10.and.iY.LE.22) THEN                                 
             dTAU(iY) = dTAU(iY-1)                                              
            ELSE                                                                
             dTAU(iY) = dTAU(2+nY-iY)                                           
           END IF                                                               
        END DO                                                                  
       ELSE                                                                     
        nY = 80                                                                 
        DO iY = 1, nY                                                           
          IF(iY.LE.11) THEN                                                     
             dTAU(iY) = 0.5*(exp(0.5*iY)-exp(0.5))                              
          ELSE IF(iY.GT.11.and.iY.LE.71) THEN                                   
             dTAU(iY) = dTAU(iY-1)                                              
          ELSE                                                                  
             dTAU(iY) = dTAU(2+nY-iY)                                           
          END IF                                                                
        END DO                                                                  
       END IF                                                                   
      END IF                                                                    

      RETURN                                                                    
      END                                                                       