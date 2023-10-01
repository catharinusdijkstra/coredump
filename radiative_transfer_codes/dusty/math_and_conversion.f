c     ======================================================================        
c     This is the file with auxiliary math and number-to-string conversion          
c     subroutines.                                             [MN, Mar,99]         
c     ======================================================================        
C         Table of Contents                                                         
C                                                                                   
C         ADD                                                                       
C         ADD2                                                                      
C         ANALINT                                                                   
C         Attach                                                                    
C         BESSEL                                                                    
C         BHMIE                                                                     
C         CHKRANGE                                                                  
C         CHKSPLIN                                                                  
C         Clean                                                                     
C         EMPTY                                                                     
C         FileMSG                                                                   
C         FindErr                                                                   
C         FINDMAX                                                                   
C         FINDRMS                                                                   
C         GETFS                                                                     
C         H                                                                         
C         KRON                                                                      
C         LINE                                                                      
C         LININTER                                                                  
C         LINSYS                                                                    
C         LUBKSB                                                                    
C         LUDCMP                                                                    
C         MAKETABLE                                                                 
C         MAPLE3                                                                    
C         MIDSQL                                                                    
C         MPROVE                                                                    
C         MSG                                                                       
C         MULTIPLY                                                                  
C         MULTIP2                                                                   
C         MYSPLINE                                                                  
C         POLINT                                                                    
C         POWERINT                                                                  
C         PRODUCT                                                                   
C         ROMBERG2                                                                  
C         ROMBY                                                                     
C         SCALETO1                                                                  
C         SHIFT                                                                     
C         SIMPSON                                                                   
C         SORT                                                                      
C         SPLINE                                                                    
C         SPLINE2                                                                   
C         SPLINT                                                                    
C         TRAPZD                                                                    
C         WriteOut                                                                  
C         ZBRAC                                                                     
C         ZRIDDR                                                                    
c     ======================================================================        
                                                                                
                                                                                
      SUBROUTINE ADD(np1,nr1,np2,nr2,q1,q2,q3,qOut)                             
c     ======================================================================        
c     This subroutine evaluates the following expression:                           
c     [qOut] = [q1] + [q2] + [q3]. qout, q1, q2 and q2 are matrices of              
c     physical size (np2,np1) and real size (nr2,nr1).     [Z.I., Nov. 1995]        
c     ======================================================================        
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER  np1, nr1, np2, nr2, i2, i1                                       
      DOUBLE PRECISION  q1(np2,np1), q2(np2,np1), q3(np2,np1),                  
     &       qOut(np2,np1)                                                      

c     loop over index 2                                                         
      DO i2 = 1, nr2                                                            
c       loop over index 1                                                       
        DO i1 = 1, nr1                                                          
          qOut(i2,i1) = q1(i2,i1) +  q2(i2,i1) + q3(i2,i1)                      
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       

                                                                                
      SUBROUTINE ADD2(flxS,flxE,fBSum,nY)                                       
c     ======================================================================        
c     This subroutine is auxiliary for finding the bolometric                       
c     components of the scattered and emitted diffuse flux.   [MN, May'99]          
c     ======================================================================        
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, iY                                                            
      DOUBLE PRECISION flxS(npL,npY), flxE(npL,npY), flxSB(npY),                
     &          flxEB(npY), fBSum(npY)                                          

      CALL Bolom(flxS,flxSB)                                                    
      CALL Bolom(flxE,flxEB)                                                    
      DO iY = 1, nY                                                             
        fBSum(iY) = flxSB(iY) + flxEB(iY)                                       
      END DO                                                                    

      RETURN                                                                    
      END                                                                       

                                                                                
      SUBROUTINE ANALINT(Nanal,xaux,yaux,m,aux,error)                           
c     ======================================================================        
c     This subroutine calculates integral I(x**m*y(x)*dx). Both y and x are         
c     1D arrays, y(i), x(i) with i=1,Nanal. The method used is approximation        
c     of y(x) by y = P(x) + d/sqrt(1-x*x), where P(x) is the polynomial of          
c     order Nanal-1, and analytic evaluation of the integral. It is assumed         
c     that xaux(1)=0. Coefficients are determined from the set of Nanal             
c     linear equations and subsequent call to the linear system solver              
c     LINSYS.                                              [Z.I., Nov. 1995]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER i, j, Nanal, error                                                
      DOUBLE PRECISION xaux(Nanal), yaux(Nanal), coeff(4), A(npY,npY),          
     &               m, aux, b                                                  

      error = 0                                                                 
c     generate matrix A and vector B                                            
      DO i = 1, Nanal                                                           
        DO j = 1, Nanal-1                                                       
          IF (xaux(i).EQ.0.0.AND.j.EQ.1) THEN                                   
            A(i,j) = 1.0                                                        
          ELSE                                                                  
            A(i,j) = xaux(i)**(1.0*j-1.0)                                       
          END IF                                                                
        END DO                                                                  
        A(i,Nanal) = 1.0/sqrt(1.0-xaux(i)*xaux(i))                              
      END DO                                                                    
                                                                                
c     solve for the coefficients                                                
      CALL LINSYS(Nanal,A,yaux,coeff,error)                                     
        IF(error.NE.0) THEN                                                     
         CALL MSG(19)                                                           
         iERROR = iERROR + 1                                                    
         RETURN                                                                 
        END IF                                                                  
c     upper limit for integration:                                              
      b = xaux(Nanal)                                                           
c     evaluate m-dependent contribution of the last term                        
      IF (m.GT.0.1) THEN                                                        
        IF (m.GT.1.1) THEN                                                      
c         this is for m=2                                                       
          aux = 0.5*(DASIN(b)-b*sqrt(1.-b*b))                                   
        ELSE                                                                    
c         this is for m=1                                                       
          aux = 1.0 - sqrt(1.-b*b)                                              
        ENDIF                                                                   
      ELSE                                                                      
c       this is for m=0                                                         
        aux = DASIN(b)                                                          
      ENDIF                                                                     
      aux = aux * coeff(Nanal)                                                  
c     add contribution from the polynom                                         
      DO i = 1, Nanal-1                                                         
        aux = aux + coeff(i) * (b**(m+1.0*i)) / (m+1.0*i)                       
      END DO                                                                    

999   RETURN                                                                    
      END                                                                       


      SUBROUTINE Attach(root,length,ext,fname)                                  
c     Attaches extensions to the root cleaned by Clean                          
c     =======================================================================       
      CHARACTER*(*) root, ext, fname                                            
      INTEGER i, length                                                         

      DO i = 1, LEN(fname)                                                      
        fname(i:i) = ' '                                                        
      END DO                                                                    
      fname(:length) = root(:length)                                            
      fname(length + 1:) = ext                                                  

      RETURN                                                                    
      END                                                                       


      SUBROUTINE Clean(StrIn, StrOut, Length)                                   
c     =======================================================================       
c     Find meaningful part of StrIn without leading and trailing junk           
c     It is returned left-justified in StrOut, right-padded with blanks         
c     The number of meaningful characters is returned in Length.                
c     In case of any problems, StrOut is empty. This sub should be used to      
c     clean every input filename immediately after DUSTY reads it. [ME,'99]     
c     =======================================================================       
      CHARACTER*(*) StrIn, StrOut                                               
      INTEGER i, first, last, Length                                            

      DO i = 1, LEN(StrOut)                                                     
        StrOut(i:i) = ' '                                                       
      END DO                                                                    
      first = 1                                                                 
      last = LEN(StrIn)                                                         
      If (first.gt.last) return                                                 
c     Find end of leading junk:                                                 
      DO WHILE (StrIn(first:first).LE.' ')                                      
       first = first + 1                                                        
       if (first.gt.last) return                                                
      END DO                                                                    
c     Find start of trailing junk:                                              
      DO WHILE (StrIn(last:last).LE.' ')                                        
       last = last - 1                                                          
       if (last.lt.first) return                                                
      END DO                                                                    
c     Now trim all junk:                                                        
      StrOut = StrIn(first:last)                                                
      Length = last - first + 1                                                 

      RETURN                                                                    
      END                                                                       


      DOUBLE PRECISION FUNCTION Bessel(x)                                       
c     =======================================================================       
c     This function evaluates the Bessel function of the zeroth kind.               
c     Formulae are from Abramowitz & Stegun.               [Z.I., Jan. 1997]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER i                                                                 
      DOUBLE PRECISION x, c(6), Pi                                              

      Pi = 2.0*ASIN(1.0)                                                        
      c(1) = -2.2499997                                                         
      c(2) =  1.2656208                                                         
      c(3) = -0.3163866                                                         
      c(4) =  0.0444479                                                         
      c(5) = -0.0039444                                                         
      c(6) =  0.00021                                                           
      Bessel=0.0                                                                
      IF (x.LE.3.0)THEN                                                         
        DO i=1,6                                                                
          Bessel = Bessel + c(i)*(x/3.0)**(2.0*i)                               
        END DO                                                                  
        Bessel = 1.0 + Bessel                                                   
        ELSE                                                                    
        Bessel = sqrt(2.0/Pi/x) * dcos(x-Pi/4.0)                                
      ENDIF                                                                     

      RETURN                                                                    
      END                                                                       


c     __________________________________________________________________                                                                             
c     This subroutine obtained from prof. P. Menguc, Dept. of Mechanical            
c     Eng., University of Kentucky.                        [Z.I., Aug. 1996]        
c
c     SUBROUTINE BHMIE CALCULATES AMPLITUDE SCATTERING MATRIX ELEMENTS          
c     & EFFICIENCIES FOR EXTINCTION, TOTAL SCATTERING AND BACSCATTERING,        
c     FOR A GIVEN SIZE PARAMETER AND RELATIVE REFRACTIVE INDEX                  
c     __________________________________________________________________        
      subroutine bhmie (x,refrel,nang,s1,s2,qext,qsca,qback)                    
      dimension amu(100),theta(100),pi(100),tau(100),pi0(100),pi1(100)          
      complex d(3000),y,refrel,xi,xi0,xi1,an,bn,s1(200),s2(200)                 
      double precision psi0,psi1,psi,dn,dx                                      
      dx=x                                                                      
      y=x*refrel                                                                
c     ___________________________________________________________________       
c     series terminated after nstop terms                                       
c     ___________________________________________________________________       
      xstop=x+4.*x**.3333 +2.0                                                  
      nstop=xstop                                                               
      ymod=cabs(y)                                                              
      nmx=amax1(xstop,ymod) + 15                                                
      dang=1.570796327/float(nang-1)                                            
      do 555 j = 1,nang                                                         
      theta(j)= (float(j)-1.)*dang                                              
555   amu(j)=cos(theta(j))                                                      

c     logarithmic derivative d(j) calculated by downward recurrence             
c     beginning with initial value 0.0 + i*0.0 at j = nmx                       
      d(nmx)=cmplx(0.0,0.0)                                                     
      nn=nmx-1                                                                  
      do 120 n=1,nn                                                             
      rn=nmx-n+1                                                                
      d(nmx-n)=(rn/y)-(1./(d(nmx-n+1)+rn/y))                                    
120   continue                                                                  
      do 666 j=1,nang                                                           
      pi0(j)=0.0                                                                
      pi1(j)=1.0                                                                
666   continue                                                                  
      nn=2*nang-1                                                               
      do 777 j=1,nn                                                             
      s1(j)=cmplx(0.0,0.0)                                                      
      s2(j)=cmplx(0.0,0.0)                                                      
777   continue                                                                  

c     riccati bessel functions with real argument x calculated by upward        
c     recurrence                                                                
      psi0=cos(dx)                                                              
      psi1=sin(dx)                                                              
      chi0=-sin(x)                                                              
      chi1=cos(x)                                                               
      apsi0=psi0                                                                
      apsi1=psi1                                                                
      xi0=cmplx(apsi0,-chi0)                                                    
      xi1=cmplx(apsi1,-chi1)                                                    
      qsca=0.0                                                                  
      n=1                                                                       
200   dn=n                                                                      
      rn=n                                                                      
      fn=(2.*rn+1.)/(rn*(rn+1.))                                                
      psi=(2.*dn-1.)*psi1/dx-psi0                                               
      apsi=psi                                                                  
      chi=(2.*rn-1.)*chi1/x -  chi0                                             
      xi = cmplx(apsi,-chi)                                                     
      an=(d(n)/refrel+rn/x)*apsi - apsi1                                        
      an=an/((d(n)/refrel+rn/x)*xi - xi1)                                       
      bn=(refrel *d(n)+rn/x)*apsi - apsi1                                       
      bn=bn/((refrel*d(n)+rn/x)*xi - xi1)                                       
      qsca=qsca+(2.*rn+1.)*(cabs(an)*cabs(an)+cabs(bn)*cabs(bn))                
      do 789 j=1,nang                                                           
      jj=2*nang-j                                                               
      pi(j)=pi1(j)                                                              
      tau(j)=rn*amu(j)*pi(j) - (rn+1.)*pi0(j)                                   
      p=(-1.)**(n-1)                                                            
      s1(j)=s1(j)+fn*(an*pi(j)+bn*tau(j))                                       
      t=(-1.)**n                                                                
      s2(j)=s2(j) + fn*(an*tau(j)+bn*pi(j))                                     
      if (j .eq. jj) go to 789                                                  
      s1(jj)=s1(jj) + fn*(an*pi(j)*p + bn*tau(j)*t)                             
      s2(jj)=s2(jj) + fn*(an*tau(j)*t + bn*pi(j)*p)                             
789   continue                                                                  
      psi0=psi1                                                                 
      psi1=psi                                                                  
      apsi1=psi1                                                                
      chi0=chi1                                                                 
      chi1=chi                                                                  
      xi1=cmplx(apsi1,-chi1)                                                    
      n=n+1                                                                     
      rn=n                                                                      
      do 999 j=1,nang                                                           
      pi1(j)=((2.*rn-1.)/(rn-1.))*amu(j)*pi(j)                                  
      pi1(j)=pi1(j) - rn*pi0(j)/(rn-1.)                                         
      pi0(j) = pi(j)                                                            
999   continue                                                                  
      if (n-1-nstop) 200, 300, 300                                              
300   qsca=(2./(x*x))*qsca                                                      
      qext=(4./(x*x))*real(s1(1))                                               
      qback=(4./(x*x))*cabs(s1(2*nang -1))*cabs(s1(2*nang -1))                  
      return                                                                    
      end                                                                       

                                                                                
      SUBROUTINE CHKRANGE(dr,x)                                                 
c     =======================================================================        
c      This subroutine checks if x is within the allowed range defined by            
c      dr<<1:                                                                        
c              dr**2 < x < 1/dr**2                                                   
c      If it is not then x = 0.0                            [Z.I., Jan. 1997]        
c     =======================================================================        
      IMPLICIT none                                                             
      DOUBLE PRECISION x, dr                                                    

      IF ((x-dr*dr)*(x-1./dr/dr).LT.0.0) THEN                                   
        continue                                                                
      ELSE                                                                      
        x = 0.0                                                                 
      END IF                                                                    

      RETURN                                                                    
      END                                                                       
                                                                                

      SUBROUTINE CHKSPLIN(x,fun,funmid,N,coef,maxerr,RDW)                       
c     ======================================================================        
c     This subroutine checks the spline coefficients coef(i,j):                     
c     fun(x)=coef(i,1) + coef(i,2)*x + coef(i,3)*x^2 + coef(i,4)*x^3,               
c     for x(i).LE.x.LE.x(i+1) with i=1..N. Array funmid(1..N-1) contains the        
c     values of function fun at mid points defined as                               
c     xmid(i)=SQRT(x(i)*x(i+1). If spline approximation produces error              
c     greater than maxerr, or funmid<0, a straight line is produced between         
c     x(i) and x(i+1).                                                              
c     ======================================================================        
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER N, i, iC                                                          
      DOUBLE PRECISION x(npY), fun(npY), funmid(npY), coef(npY,4),              
     &       maxerr, error, slope, xmid, funSpline, aux, power, yR, yL          
      LOGICAL RDW                                                               

c       check the midpoints                                                     
        DO i = 1, N - 1                                                         
          xmid = dsqrt(x(i)*x(i+1))                                             
          funSpline = 0.0                                                       
          DO iC=1,4                                                             
            IF (xmid.EQ.0.0.AND.iC.EQ.1) THEN                                   
              aux = 1.0                                                         
              ELSE                                                              
              aux = xmid**(float(iC)-1.0)                                       
            END IF                                                              
            funSpline = funSpline + coef(i,iC)*aux                              
          END DO                                                                
          error = DABS((funSpline-funmid(i))/funmid(i))                         
c         check for the deviation at the midpoint                               
          IF (error.GE.maxerr.OR.funSpline.LE.0.0) THEN                         
            slope = (fun(i+1) - fun(i)) / (x(i+1)-x(i))                         
            coef(i,1) = fun(i) - x(i) * slope                                   
            coef(i,2) = slope                                                   
            coef(i,3) = 0.0                                                     
            coef(i,4) = 0.0                                                     
          END IF                                                                
c         check for the logarithmic derivative (only for RDW)                   
          IF(RDW) THEN                                                          
            yL = fun(i)                                                         
            yR = fun(i+1)                                                       
            IF (x(i)*x(i+1).GT.0.AND.yL*yR.GT.0) THEN                           
              power = log(yR/yL)/log(x(i+1)/x(i))                               
              IF (abs(power).GT.10.) THEN                                       
                slope = (yR - yL) / (x(i+1)-x(i))                               
                coef(i,1) = yL - x(i) * slope                                   
                coef(i,2) = slope                                               
                coef(i,3) = 0.0                                                 
                coef(i,4) = 0.0                                                 
              END IF                                                            
            END IF                                                              
          END IF                                                                
        END DO                                                                  

      RETURN                                                                    
      END                                                                       


      SUBROUTINE CHKSPLINold(x,fun,funmid,N,coef,maxerr)                        
c     =======================================================================       
c     This subroutine checks the spline coefficients coef(i,j):                     
c     fun(x)=coef(i,1) + coef(i,2)*x + coef(i,3)*x^2 + coef(i,4)*x^3,               
c     for x(i).LE.x.LE.x(i+1) with i=1..N. Array funmid(1..N-1) contains the        
c     values of function fun at mid points defined as                               
c     xmid(i)=SQRT(x(i)*x(i+1). If spline approximation produces error              
c     greater than maxerr, or funmid<0, a straight line is produced between         
c     x(i) and x(i+1).                                     [Z.I., Feb. 1995]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER N, i, iC                                                          
      DOUBLE PRECISION x(npY), fun(npY), funmid(npY), coef(npY,4), aux,         
     &     funSpline, maxerr, error, slope, xmid                                

c       check the midpoints                                                     
        DO i = 1, N - 1                                                         
          xmid = dsqrt(x(i)*x(i+1))                                             
          funSpline = 0.0                                                       
          DO iC=1,4                                                             
            IF (xmid.EQ.0.0.AND.iC.EQ.1) THEN                                   
              aux = 1.0                                                         
              ELSE                                                              
              aux = xmid**(float(iC)-1.0)                                       
            END IF                                                              
            funSpline = funSpline + coef(i,iC)*aux                              
          END DO                                                                
          error = DABS((funSpline-funmid(i))/funmid(i))                         
          IF (error.GE.maxerr.OR.funSpline.LE.0.0) THEN                         
            slope = (fun(i+1) - fun(i)) / (x(i+1)-x(i))                         
            coef(i,1) = fun(i) - x(i) * slope                                   
            coef(i,2) = slope                                                   
            coef(i,3) = 0.0                                                     
            coef(i,4) = 0.0                                                     
          END IF                                                                
        END DO                                                                  

      RETURN                                                                    
      END                                                                       


      INTEGER FUNCTION EMPTY(line)                                              
c     =======================================================================       
c     This function is 1 if string 'line' is empty, or if it contains only          
c     '%', and 0 otherwise.                                                         
c                                                          [Z.I., Nov. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER i, iTeX, l                                                        
      CHARACTER ch                                                              
      CHARACTER*(*) line                                                        

      l = LEN(line)                                                             
      EMPTY = 1                                                                 
      iTeX = 0                                                                  
      DO i = 1, l                                                               
        ch = line(i:i)                                                          
        IF(EMPTY.EQ.1.AND.ch.EQ.'%') iTeX = 1                                   
         IF (ch.NE.' ') EMPTY = 0                                               
      END DO                                                                    
      IF (iTeX.EQ.1) EMPTY = 1                                                  

      RETURN                                                                    
      END                                                                       


      SUBROUTINE FileMSG(fname,strg)                                            
c     =======================================================================       
c         Prints a message in *.out file in case of error opening the user          
c         supplied files.                                                           
c     =======================================================================       
      IMPLICIT NONE                                                             
      CHARACTER aux*230, strg*(*), fname*(*)                                    
      INTEGER length, Empty                                                     

1     read(1,'(a)') aux                                                         
      IF (Empty(aux).EQ.1) goto 1                                               
      CALL Clean(aux,fname,length)                                              
                                                                                
      open(10, ERR=100, FILE=fname, STATUS='OLD') 	                             
      close(10)                                                                 
	RETURN                                                                         
                                                                                
100   write(12,*)' *** FATAL ERROR IN DUSTY! **************************'        
      write(12,*)' File with the ',strg                                         
      write(12,'(a2,a)')'  ',fname                                              
      write(12,*)' is missing ?!'                                               
      write(12,*)' ****************************************************'        
      close(12)                                                                 

      STOP                                                                      
      END                                                                       


      SUBROUTINE FindErr(fbol,maxFerr,nY)                                       
c     ========================================================================      
c     This subroutine finds maximum err in flux conservation for both               
c     spherical and slab case as (fmax-fmin)/(fmax+fmin)   [MN,Aug'99]              
c     =========================================================================     
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER nY, iY                                                            
      DOUBLE PRECISION fbol(npY), maxFerr, fmin, fmax, aux                      

c     Find the min and max of fbol values                                       
c     The abs and lower limit on fbol are protection for the case               
c     of completely symmetric slab illumination. The lower limit                
c     is bound by the numerical accuracy of the flux calculation                
        fmin = 1.e5                                                             
        fmax = 0.                                                               
        DO iY = 1, nY                                                           
           aux = fbol(iY)                                                       
           IF (ksi.eq.1.0) aux = dabs(aux)                                      
           IF (dabs(aux).LE.accFbol) aux = accFbol                              
           IF(aux.LT.fmin) fmin = aux                                           
           IF(aux.GT.fmax) fmax = aux                                           
        END DO                                                                  
        if (fmax.LT.0.) then                                                    
c     bad solution; overall flux cannot be negative                             
            maxFerr = 1                                                         
        else                                                                    
            maxFerr = (fmax - fmin)/(fmax + dabs(fmin))                         
        end if                                                                  

      RETURN                                                                    
      END                                                                       
                                                                                

      SUBROUTINE FindMax(NN,i1,i2,A,Amax)                                       
c     =======================================================================       
c     This subroutine finds maximum values, Amax, of an array A(nY) between         
c     values A(i1) and A(i2).                              [Z.I., Jul. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER NN, i1, i2, i                                                     
      DOUBLE PRECISION A(NN), Amax                                              

      Amax = A(i1)                                                              
c     loop over radial positions                                                
      DO i = i1, i2                                                             
        IF (A(i).GT.Amax) Amax = A(i)                                           
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE FindRMS(typ,X,val,accur,N)                                     
c     =======================================================================       
c     Finds relative deviations 'accur' of an array X(N) from a given value val.    
c     For typ=1 accur is maximal deviation, and for typ=2 the rms deviation.        
c                                                             [ZI'95; MN'99]        
c     =======================================================================       
      IMPLICIT NONE                                                             
      INTEGER N, i, typ                                                         
      DOUBLE PRECISION X(N), val, accur, ss, dev                                

      IF (typ.EQ.1) THEN                                                        
        accur = 0.0                                                             
        DO i = 1, N                                                             
          dev = (X(i)-val)/val                                                  
          IF (DABS(dev).GT.accur) accur = DABS(dev)                             
        END DO                                                                  
      ELSE                                                                      
        ss = 0.0                                                                
        DO i = 1, N                                                             
          dev = X(i)-val                                                        
          ss = ss + dev*dev                                                     
        END DO                                                                  
        accur = sqrt(ss/N/(N-1.))                                               
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE GetFS(xx,nm,flag,str)                                          
c     =======================================================================       
c     This subroutine writes number xx to a string str according to a format        
c     f?.nm. Here ? stands for the number of needed places. A blank is              
c     inserted at the beginning, and for flag.NE.1 another one if number is         
c     positive. If xx<0 second blank is replaced by '-'. For example, for           
c     flag=0 and xx = -0.1234E+02, calling this subroutine with nm=1 will           
c     result in str = ' -12.3', while xx = 0.0123 with nm=3 gives '  0.012'.        
c     If flag=1 minus will be ignored, for example xx = -0.1234E+02 and nm=1        
c     will result in str = ' 12.3',                        [Z.I., Nov. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      CHARACTER ch                                                              
      CHARACTER*(*) str                                                         
      INTEGER  flag, nm, db, i, d(20), j, k, dnmp1                              
      DOUBLE PRECISION xx, x, rest                                              

      DO i = 1, len(str)                                                        
        str(i:i) = ' '                                                          
      END DO                                                                    
                                                                                
      x = xx                                                                    
      str(1:1) = ' '                                                            
      i = 2                                                                     
      IF (flag.NE.1) THEN                                                       
         IF (x.LT.0.0) THEN                                                     
           str(i:i) = '-'                                                       
         ELSE                                                                   
           str(i:i) = ' '                                                       
         END IF                                                                 
         i = i + 1                                                              
      END IF                                                                    
      IF (x.LT.0.0) x = -x                                                      
c     first check if x will have to be rounded up                               
c     find (nm+1)-th decimal digit                                              
      dnmp1 = int(x*10.**(nm+1)-int(x*10.**nm)*10.)                             
      IF (dnmp1.GE.5) x = x + 1./10.0**nm                                       
      IF (x.GE.1.0) THEN                                                        
c       number of digits before the decimal sign                                
        db = int(log10(x)) + 1                                                  
c       copy all these digits to str                                            
        DO j = 1, db                                                            
          rest = x                                                              
          IF (j.GT.1) THEN                                                      
            DO k = 1, j-1                                                       
              rest = rest - d(k)*10.**(db-k)                                    
            END DO                                                              
          END IF                                                                
          d(j) = int(rest/10.**(db-j))                                          
          write(ch,'(i1)')d(j)                                                  
          str(i:i) = ch                                                         
          i = i + 1                                                             
        END DO                                                                  
        rest = rest - d(db)                                                     
        IF (nm.GT.0) THEN                                                       
          str(i:i) = '.'                                                        
          i = i + 1                                                             
        END IF                                                                  
      ELSE                                                                      
        str(i:i) = '0'                                                          
        i = i + 1                                                               
        IF (nm.GT.0) THEN                                                       
          str(i:i) = '.'                                                        
          i = i + 1                                                             
        END IF                                                                  
        rest = x                                                                
      END IF                                                                    
c     now copy all nm remaining decimal digits to str                           
      IF (nm.GT.0) THEN                                                         
        DO j = 1, nm                                                            
          d(j) = int(rest*10.**j)                                               
          IF (j.GT.1) THEN                                                      
            DO k = 1, j-1                                                       
              d(j)=d(j)-int(d(k)*10.**(j-k))                                    
            END DO                                                              
          END IF                                                                
          write(ch,'(i1)')d(j)                                                  
          str(i:i) = ch                                                         
          i = i + 1                                                             
        END DO                                                                  
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      DOUBLE PRECISION FUNCTION H(x1,x2)                                        
c     =======================================================================       
c     This function calculates the step function: H=1 for x1 >= x2 and H=0          
c     for x1 < x2.                                         [Z.I., Nov. 1995]        
c     =======================================================================       
      IMPLICIT none                                                             
      DOUBLE PRECISION x1, x2                                                   

      IF (x1.GE.x2) THEN                                                        
        H = 1.0                                                                 
      ELSE                                                                      
        H = 0.0                                                                 
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      INTEGER FUNCTION Kron(i1,i2)                                              
c     =======================================================================       
c     This function is Kronecker delta-function defined as:                         
c     Kron(i1,i2) = 1 for i1=i2                                                     
c     Kron(i1,i2) = 0 otherwise.                           [Z.I., Dec. 1995]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER i1, i2                                                            

      IF (i1.EQ.i2) THEN                                                        
        Kron = 1                                                                
      ELSE                                                                      
        Kron = 0                                                                
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE LINE(com,typ,unt)                                              
c     =======================================================================       
c     This subroutine writes a line into file open as unt. For type = 1             
c     the line is '---', and for type = 2 '==='.If com=1 a comment sign # is        
c     added in the beginning (this is when line is used in file headers)            
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER com, typ, unt                                                     

      IF(typ.EQ.1) THEN                                                         
       IF(com.eq.1) THEN                                                        
        write(unt,'(a50)')                                                      
     &   '# ------------------------------------------------'                   
       ELSE                                                                     
        write(unt,*)'--------------------------------------------------'        
       END IF                                                                   
      ELSE                                                                      
       IF(com.eq.1) THEN                                                        
        write(unt,'(a50)')                                                      
     &   '# ================================================'                   
       ELSE                                                                     
        write(unt,*)'=================================================='        
       END IF                                                                   
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE LinInter(NN,N,x,y,xloc,iNloc,yloc)                             
c     =======================================================================       
c     This subroutine performs linear interpolation for y(x) such that              
c     yloc = y(xloc). It is assumed that x is monotonously increasing.              
c                                                          [Z.I., Mar. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER NN, N, i, istop, iNloc                                            
      DOUBLE PRECISION x(NN), y(NN), xloc, yloc                                 

      IF (N.GT.1) THEN                                                          
        IF ((x(1)-xloc)*(x(N)-xloc).LE.0.0) THEN                                
          istop = 0                                                             
          i = 1                                                                 
          DO WHILE (istop.NE.1)                                                 
            i = i + 1                                                           
            IF (i.GT.N) stop 'LinInter ???'                                     
            IF (x(i).GE.xloc) THEN                                              
              istop = 1                                                         
              iNloc = i                                                         
              yloc = y(i-1) + (y(i)-y(i-1))/(x(i)-x(i-1))*(xloc-x(i-1))         
            END IF                                                              
          END DO                                                                
          ELSE                                                                  
          IF (xloc.LE.x(1)) yloc = y(1)                                         
          IF (xloc.GE.x(N)) yloc = y(N)                                         
        END IF                                                                  
        ELSE                                                                    
        yloc = y(1)                                                             
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE LINSYS(Nreal,A,B,X,error)                                      
c     =======================================================================       
c     This subroutine solves the set of linear equations [A]*[X] = [B] for          
c     X [A(k,1)*X(1)+A(k,2)*X(2)+...+A(k,Nreal)*X(Nreal) = B(k), k=1,Nreal).        
c     The real size of matrix A is Nreal x Nreal and its physical dimension         
c     is npY x npY, where npY comes from INCLUDE 'userpar.inc'. Both vectors        
c     B and X have real lengths Nreal. The set is solved by calls to LUDCMP         
c     and LUBKSB and the solution is improved subsequently by a call to             
c     MPROVE. These three subroutines are taken from Numerical Recipes.             
c                                                          [Z.I., Nov. 1995]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER Nreal, indx(npY), i, j, error                                     
      DOUBLE PRECISION A(npY,npY), B(npY), X(npY)                               
      DOUBLE PRECISION A1(npY,npY), B1(npY), A2(npY,npY), B2(npY), d            

      error = 0                                                                 
c     generate DOUBLE PRECISION copies of A and B (two copies because they          
c     are changed in LUDCMP and LUBKSB, but still needed for MPROVE)                
      DO i = 1, Nreal                                                           
        B1(i) = B(i)                                                            
        B2(i) = B(i)                                                            
        DO j = 1, Nreal                                                         
           A1(i,j) = A(i,j)                                                     
           A2(i,j) = A(i,j)                                                     
        END DO                                                                  
      END DO                                                                    
c     solve the system                                                          
      CALL LUDCMP(A1,Nreal,npY,indx,d,error)                                    
      IF (error.NE.0) RETURN                                                    
      CALL LUBKSB(A1,Nreal,npY,indx,B1)                                         
c     improve the solution (saved in B)                                         
      CALL MPROVE(A2,A1,Nreal,npY,indx,B2,B1)                                   
c     copy the improved solution to output vector X                             
      DO i = 1, Nreal                                                           
        X(i) = B1(i)                                                            
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE LUBKSB(A,N,NP,INDX,B)                                          
c     =======================================================================       
      DIMENSION INDX(NP)                                                        
      DOUBLE PRECISION A(NP,NP),B(NP)                                           

      II=0                                                                      
      DO 12 I=1,N                                                               
      LL=INDX(I)                                                                
      SUM=B(LL)                                                                 
      B(LL)=B(I)                                                                
      IF (II.NE.0)THEN                                                          
        DO 11 J=II,I-1                                                          
          SUM=SUM-A(I,J)*B(J)                                                   
11        CONTINUE                                                              
      ELSE IF (SUM.NE.0.) THEN                                                  
        II=I                                                                    
      ENDIF                                                                     
      B(I)=SUM                                                                  
12    CONTINUE                                                                  
      DO 14 I=N,1,-1                                                            
      SUM=B(I)                                                                  
      IF(I.LT.N)THEN                                                            
        DO 13 J=I+1,N                                                           
          SUM=SUM-A(I,J)*B(J)                                                   
13        CONTINUE                                                              
      ENDIF                                                                     
      B(I)=SUM/A(I,I)                                                           
14    CONTINUE                                                                  

      RETURN                                                                    
      END                                                                       
                                                                                

      SUBROUTINE LUDCMP(A,N,NP,INDX,D,error)                                    
c     =======================================================================       
      PARAMETER (NMAX=10000,TINY=1.0E-20)                                       
      DIMENSION INDX(NP)                                                        
      INTEGER error                                                             
      DOUBLE PRECISION A(NP,NP),VV(NMAX), D, SUM                                

      error = 0                                                                 
      D = 1.                                                                    
      DO I = 1, N                                                               
       AAMAX=0.                                                                 
       DO J = 1, N                                                              
        IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))                           
       END DO                                                                   
       IF (AAMAX.EQ.0.) THEN                                                    
        error = 5                                                               
        RETURN                                                                  
       ENDIF                                                                    
       VV(I)=1./AAMAX                                                           
      END DO                                                                    
      DO J = 1 , N                                                              
       IF (J.GT.1) THEN                                                         
        DO I = 1, J-1                                                           
          SUM=A(I,J)                                                            
          IF (I.GT.1)THEN                                                       
            DO K = 1, I-1                                                       
             SUM=SUM-A(I,K)*A(K,J)                                              
            END DO                                                              
            A(I,J)=SUM                                                          
          ENDIF                                                                 
        END DO                                                                  
       ENDIF                                                                    
       AAMAX=0.                                                                 
       DO I = J, N                                                              
        SUM=A(I,J)                                                              
        IF (J.GT.1)THEN                                                         
          DO K = 1, J-1                                                         
            SUM=SUM-A(I,K)*A(K,J)                                               
          END DO                                                                
          A(I,J)=SUM                                                            
        ENDIF                                                                   
        DUM=VV(I)*DABS(SUM)                                                     
        IF (DUM.GE.AAMAX) THEN                                                  
          IMAX=I                                                                
          AAMAX=DUM                                                             
        ENDIF                                                                   
       END DO                                                                   
       IF (J.NE.IMAX)THEN                                                       
        DO K = 1, N                                                             
          DUM=A(IMAX,K)                                                         
          A(IMAX,K)=A(J,K)                                                      
          A(J,K)=DUM                                                            
        END DO                                                                  
        D=-D                                                                    
        VV(IMAX)=VV(J)                                                          
       ENDIF                                                                    
       INDX(J)=IMAX                                                             
       IF(J.NE.N)THEN                                                           
        IF(A(J,J).EQ.0.)A(J,J)=TINY                                             
        DUM=1./A(J,J)                                                           
        DO I = J+1, N                                                           
          A(I,J)=A(I,J)*DUM                                                     
        END DO                                                                  
       ENDIF                                                                    
      END DO                                                                    
      IF(A(N,N).EQ.0.)A(N,N)=TINY                                               

      RETURN                                                                    
      END                                                                       


      SUBROUTINE MakeTable(Elems,rows,cols,unt)                                 
c     =======================================================================       
c         This is an auxiliary subroutine for print out of tables                   
c         of Elems(cols,rows) in output unit 'unt'. rows = max{npL,npY}.            
c         This array is defined in PrOut as well and if you change its size         
c         make sure you do this in both places.                 [MN, Mar'98]        
c     =======================================================================       
      IMPLICIT NONE                                                             
      INTEGER rows, cols, unt, k, i                                             
      DOUBLE PRECISION Elems(25,200)                                            

      DO i = 1, rows                                                            
        write(unt,'(1p,21E11.3)') (Elems(k,i),k=1,cols)                         
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE Maple3(w,z,p,MpInt)                                            
c     =======================================================================       
c     This function calculates indefinite integral:                                 
c        MpInt(iC) = INT(w^(2-iC) / sqrt(w^2-p^2) * dw), for iC=1,2,3,4.            
c                                                          [Z.I., Apr. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      DOUBLE PRECISION w, z, p, MpInt(4)                                        

c     integrals                                                                 
      MpInt(1) = z                                                              
      MpInt(2) = dlog(w+z)                                                      
      IF (p.GT.0.0) THEN                                                        
        MpInt(3) = dacos(p/w)/p                                                 
        MpInt(4) = z/w/p/p                                                      
        ELSE                                                                    
        MpInt(3) = -1.0 / w                                                     
        MpInt(4) = -0.5 / w / w                                                 
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE midsql(funk,aa,bb,s,n)                                         
c     =======================================================================       
      INTEGER n                                                                 
      DOUBLE PRECISION aa,bb,s,funk                                             
      EXTERNAL funk                                                             
      INTEGER it,j                                                              
      DOUBLE PRECISION ddel,del,sum,tnm,x,func,a,b                              

      func(x)=2.*x*funk(aa+x**2)                                                
      b=dsqrt(bb-aa)                                                            
      a=0.                                                                      
      if (n.eq.1) then                                                          
        s=(b-a)*func(0.5*(a+b))                                                 
      else                                                                      
        it=3**(n-2)                                                             
        tnm=it                                                                  
        del=(b-a)/(3.*tnm)                                                      
        ddel=del+del                                                            
        x=a+0.5*del                                                             
        sum=0.                                                                  
        do 11 j=1,it                                                            
          sum=sum+func(x)                                                       
          x=x+ddel                                                              
          sum=sum+func(x)                                                       
          x=x+del                                                               
11      continue                                                                
        s=(s+(b-a)*sum/tnm)/3.                                                  
      endif                                                                     

      RETURN                                                                    
      END                                                                       


      SUBROUTINE MPROVE(A,ALUD,N,NP,INDX,B,X)                                   
c     =======================================================================       
      PARAMETER (NMAX=10000)                                                    
      DIMENSION INDX(N)                                                         
      DOUBLE PRECISION SDP,A(NP,NP),ALUD(NP,NP),B(N),X(N),R(NMAX)               

      DO i = 1, N                                                               
       SDP = -B(i)                                                              
       DO j = 1, N                                                              
         SDP = SDP + A(i,j)*X(j)                                                
       END DO                                                                   
       R(i) = SDP                                                               
      END DO                                                                    
      CALL LUBKSB(ALUD,N,NP,INDX,R)                                             
      DO i = 1, N                                                               
       X(i) = X(i) - R(i)                                                       
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE MSG(msgno)                                                     
c     =======================================================================       
c     This subroutine writes runtime messages to auxiliary file fname.m##           
c     or to the output file fname.out.             [ZI,Feb'96; MN,Jul'99]           
c     =======================================================================       
      IMPLICIT none                                                             
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER  msgno                                                            

      IF (msgno.EQ.1.AND.iX.GT.0) THEN                                          
       write(18,*)' ************  WARNING  ************'                        
       write(18,*)' Temperature calculation in FindTemp'                        
       write(18,*)' achieved the limit of 500 iterations'                       
      END IF                                                                    
      IF (msgno.EQ.2.AND.iX.GT.0) THEN                                          
       write(18,*)' ************  WARNING  ************'                        
       write(18,*)' Energy density iterations in RADTRANSF'                     
       write(18,*)' achieved the limit of 10000 iterations'                     
      END IF                                                                    
      IF (msgno.EQ.3) THEN                                                      
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' * Denstyp is not between 1 and 7!*'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.4.AND.iX.GT.0) THEN                                          
       write(18,*)' ************  WARNING  *****************'                   
       write(18,*)' Could not bracket in Zbrac (in sub FindTemp)'               
       write(18,*)' Something might be wrong in your input.     '               
      END IF                                                                    
      IF (msgno.EQ.5) THEN                                                      
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' *  All abundances must be >= 0!  *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.6) THEN                                                      
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' * Wavelengths for the power-law  *'                         
       write(12,*)' * spectrum must be ascending!    *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.7) THEN                                                      
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' * Relative luminosities must add *'                         
       write(12,*)' * up to a number >0!!!           *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.8) THEN                                                      
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' *    A black body temperature    *'                         
       write(12,*)' *        should be > 0 !!!       *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.9) THEN                                                      
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' * Flag for optical properties    *'                         
       write(12,*)' * should be between 1 and 3!!!   *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.10) THEN                                                     
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' * Flag for size distribution     *'                         
       write(12,*)' * should be between 1 and 3!!!   *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.11) THEN                                                     
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' * The flag for external spectrum *'                         
       write(12,*)' * should be between 1 and 6 !!!  *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.12) THEN                                                     
       write(12,*)' ***  FATAL ERROR IN DUSTY  **********'                      
       write(12,*)' Only three types of the point spread '                      
       write(12,*)' function are allowed: 1, 2 or 3 !!!  '                      
       write(12,*)' Check input file and try again       '                      
       write(12,*)' ***********************************'                        
      END IF                                                                    
c     msg 14 is not called in this version.                                     
      IF (msgno.EQ.14.AND.iX.GT.0) THEN                                         
       write(18,*)' ******** MESSAGE FROM SLBSolve *******'                     
       write(18,*)' Convergence on en.density is too slow.'                     
       write(18,*)' If the accuracy is not reached yet    '                     
       write(18,*)' will increase grid size and try again.'                     
       write(18,*)' **************************************'                     
      END IF                                                                    
      IF (msgno.EQ.15) THEN                                                     
       write(12,*) ' **************** WARNING ******************'               
       write(12,*) '  NO calculation in this case. Parameter npY'               
       write(12,*) '  needs to be at least 50. Use of the slab  '               
       write(12,*) '  parameters is suggested (see userpar.inc) '               
       write(12,*) ' *******************************************'               
      END IF                                                                    
      IF (msgno.EQ.16) THEN                                                     
       write(12,*)' ****************  WARNING  ******************'              
       write(12,*)'  The density profile Eta is too steep and the'              
       write(12,*)'  code can not handle this. Try decreasing the'              
       write(12,*)'  outer radius Y (see Manual, 3.3.3).         '              
       write(12,*)' *********************************************'              
       IF(iX.GT.0) THEN                                                         
       write(18,*)' ****************  WARNING  ******************'              
       write(18,*)'  The density profile Eta is too steep and the'              
       write(18,*)'  code can not handle this. Try decreasing the'              
       write(18,*)'  outer radius Y (see Manual, 3.3.3).         '              
       write(18,*)' *********************************************'              
       END IF                                                                   
      END IF                                                                    
      IF (msgno.EQ.17) THEN                                                     
       write(12,*)' *****************  WARNING  ********************'           
       write(12,*)'  Eta is too steep and reaches values less than  '           
       write(12,*)'  1e-12. Try decreasing the outer radius Y.      '           
       write(12,*)'  (see Manual,3.3.3)                             '           
       write(12,*)' ************************************************'           
       IF(iX.GT.0) THEN                                                         
       write(18,*)' *****************  WARNING  ********************'           
       write(18,*)'  Eta is too steep and reaches values less than  '           
       write(18,*)'  1e-12. Try decreasing the outer radius Y.      '           
       write(18,*)'  (see Manual,3.3.3)                             '           
       write(18,*)' ************************************************'           
       END IF                                                                   
      END IF                                                                    
      IF (msgno.EQ.18) THEN                                                     
       write(12,*)' ************  WARNING  ************************ '           
       write(12,*)'  Dynamical range of Eta more than 1.E-12.       '           
       write(12,*)'  The outer radius Y must be decreased so that   '           
       write(12,*)'  Eta does not go below 1.E-12 (see Manual,3.3.3)'           
       write(12,*)' *********************************************** '           
       IF(iX.GT.0) THEN                                                         
       write(18,*)' ************  WARNING  ************************ '           
       write(18,*)'  Dynamical range of Eta more than 1.E-12.       '           
       write(18,*)'  The outer radius Y must be decreased so that   '           
       write(18,*)'  Eta does not go below 1.E-12 (see Manual,3.3.3)'           
       write(18,*)' *********************************************** '           
       END IF                                                                   
      END IF                                                                    
      IF (msgno.EQ.19) THEN                                                     
       write(12,*)' ************ A BIG ERROR!!!************* '                  
       write(12,*)'  Singular matrix in LUDCMP when called   '                  
       write(12,*)'  from ANALINT. Stopping the calculation. '                  
       write(12,*)' **************************************** '                  
      END IF                                                                    
      IF (msgno.EQ.20) THEN                                                     
       write(12,*)' ************ A BIG ERROR!!! ************ '                  
       write(12,*)'  Singular matrix in LUDCMP when called   '                  
       write(12,*)'  from INVERT. Stopping the calculation.  '                  
       write(12,*)' **************************************** '                  
      END IF                                                                    

101   RETURN                                                                    
      END                                                                       


      SUBROUTINE MULTIPLY(type,np1,nr1,np2,nr2,mat,vec1,omat,flag,q1,q2)        
c     =======================================================================       
c     This subroutine evaluates the following expression:                           
c     [q2] = flag*[q1] + [mat]*[tt*vec1]. Here tt is [omat] for type=1 and          
c     1-[omat] for type=2. mat is matrix of physical size (np2,np1,np1) and         
c     real size (nr2,nr1,nr1). omat, vec1, q1 and q2 are matrices of                
c     physical size (np2,np1) and real size (nr2,nr1).     [Z.I., Nov. 1995]        
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
      INTEGER type, np1, nr1, np2, nr2, flag, i2, i1, idum                      
      DOUBLE PRECISION mat(np2,np1,np1), vec1(np2,np1), omat(np2,np1),          
     &       aux, q1(np2,np1), q2(np2,np1)                                      

c     loop over index 2                                                         
      DO i2 = 1, nr2                                                            
c       loop over index 1                                                       
        DO i1 = 1, nr1                                                          
          q2(i2,i1) = flag * q1(i2,i1)                                          
c         loop over dummy index (multiplication)                                
          DO idum = 1, nr1                                                      
            IF (type.EQ.1) THEN                                                 
              aux = omat(i2,idum)                                               
              ELSE                                                              
              aux = 1.0 - omat(i2,idum)                                         
            END IF                                                              
            q2(i2,i1) = q2(i2,i1) + mat(i2,i1,idum)*aux*vec1(i2,idum)           
          END DO                                                                
          IF (q2(i2,i1).LT.dynrange*dynrange) q2(i2,i1) = 0.0                   
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE MULTIP2(type,np1,nr1,np2,nr2,nr3,np3,mat,vec1,omat,q1)         
c     =======================================================================       
c     This subroutine evaluates the following expression:                           
c     [q1] = [mat]*[tt*vec1] / 4Pi. Here tt is [omat] for type=1 and                
c     1-[omat] for type=2. mat is matrix of physical size (np2,np3,np1) and         
c     real size (nr2,nr3,nr1). omat and vec1 are matrices of physical size          
c     (np2,np1) and real size (nr2,nr1). q1 is a matrix of physical size            
c     (np2,np3) and real size (nr2,nr3)                                             
c     1, 2 and 3 correspond to nY, nL and nP.              [Z.I., Nov. 1995]        
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
      INTEGER type, np1, nr1, np2, nr2, np3, nr3, i2, i3, idum                  
      DOUBLE PRECISION mat(np2,np3,np1), vec1(np2,np1), omat(np2,np1),          
     &       aux, q1(np2,np3)                                                   

c     loop over index 2 (wavelength)                                            
      DO i2 = 1, nr2                                                            
c       loop over index 3 (impact parameter)                                    
        DO i3 = 1, nr3                                                          
          q1(i2,i3) = 0.0                                                       
c         loop over dummy index (multiplication)                                
          DO idum = 1, nr1                                                      
            IF (type.EQ.1) THEN                                                 
              aux = omat(i2,idum)                                               
            ELSE                                                                
              aux = 1.0 - omat(i2,idum)                                         
            END IF                                                              
            q1(i2,i3) = q1(i2,i3) + mat(i2,i3,idum)*aux*vec1(i2,idum)           
          END DO                                                                
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE MYSPLINE(x,N,alpha,beta,gamma,delta)                           
c     =======================================================================       
c     This subroutine finds arrays alpha, beta, gamma and delta describing          
c     a cubic spline approximation of an unknown function f(x) given as an          
c     array f(i)=f(x(i)) with i=1..N. The cubic spline approximation is:            
c     f(x)=a(i) + b(i)*t + c(i)*t^2 + d(i)*t^3  for x(i).LE.x.LE.x(i+1)             
c     and t = (x-x(i))/(x(i+1)-x(i)), i=1..N-1. Coefficients a,b,c,d are            
c     equal to:                                                                     
c     a(i) = alpha(i,1)*f(1) + alpha(i,2)*f(2) + ... + alpha(i,N)*f(N)              
c     and b,c,d analogously.                               [Z.I., Dec. 1995]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER N, i, j, dummy, Kron                                              
      DOUBLE PRECISION x(npY), alpha(npY,npY), beta(npY,npY),                   
     &       delta(npY,npY), secnder(npY,npY), yaux(npY), deraux(npY),          
     &       y2at1, y2atN, D, gamma(npY,npY)                                    
      EXTERNAL Kron                                                             

c     generate second derivatives, secnder(j,l)                                 
      DO j = 1, N                                                               
        DO dummy = 1, N                                                         
          IF (dummy.EQ.j) THEN                                                  
            yaux(dummy) = 1.0                                                   
            ELSE                                                                
            yaux(dummy) = 0.0                                                   
          END IF                                                                
        END DO                                                                  
        y2at1 = (yaux(2)-yaux(1))/(x(2)-x(1))                                   
        y2atN = (yaux(N)-yaux(N-1))/(x(N)-x(N-1))                               
        CALL SPLINE(x,yaux,N,y2at1,y2atN,deraux)                                
        DO i = 1, N                                                             
          secnder(i,j) =  deraux(i)                                             
c          secnder(i,j) = 0.0                                                   
        END DO                                                                  
      END DO                                                                    
c     generate alpha, beta, gamma, delta                                        
      DO i = 1, N-1                                                             
        D = (x(i+1) - x(i))*(x(i+1) - x(i)) / 6.0                               
        DO j = 1, N                                                             
          alpha(i,j) = Kron(i,j)*1.0                                            
          beta(i,j) = Kron(i+1,j) - Kron(i,j)                                   
          beta(i,j) = beta(i,j) - D*(2.*secnder(i,j)+secnder(i+1,j))            
          gamma(i,j) = 3. * D * secnder(i,j)                                    
          delta(i,j) = D*(secnder(i+1,j)-secnder(i,j))                          
        END DO                                                                  
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE polint(xa,ya,n,x,y,dy)                                         
c     =======================================================================       
      INTEGER n,NMAX                                                            
      DOUBLE PRECISION dy,x,y,xa(n),ya(n)                                       
      PARAMETER (NMAX=10)                                                       
      INTEGER i,m,ns                                                            
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)                     

      ns=1                                                                      
      dif=DABS(x-xa(1))                                                         
      do 11 i=1,n                                                               
        dift=DABS(x-xa(i))                                                      
        if (dift.lt.dif) then                                                   
          ns=i                                                                  
          dif=dift                                                              
        endif                                                                   
        c(i)=ya(i)                                                              
        d(i)=ya(i)                                                              
11    continue                                                                  
      y=ya(ns)                                                                  
      ns=ns-1                                                                   
      do 13 m=1,n-1                                                             
        do 12 i=1,n-m                                                           
          ho=xa(i)-x                                                            
          hp=xa(i+m)-x                                                          
          w=c(i+1)-d(i)                                                         
          den=ho-hp                                                             
          if(den.eq.0.)pause 'failure in polint'                                
          den=w/den                                                             
          d(i)=hp*den                                                           
          c(i)=ho*den                                                           
12      continue                                                                
        if (2*ns.lt.n-m)then                                                    
          dy=c(ns+1)                                                            
        else                                                                    
          dy=d(ns)                                                              
          ns=ns-1                                                               
        endif                                                                   
        y=y+dy                                                                  
13    continue                                                                  

      RETURN                                                                    
      END                                                                       
                                                                                

      SUBROUTINE PowerInt(N,N1,N2,x,y,integral)                                 
c     =======================================================================       
c     This subroutine calculates integral I(y(x)*dx). Both y and x are              
c     1D arrays, y(i), x(i) with i=1,N (declared with NN). Lower and upper          
c     integration limits are x(N1) and x(N2), respectively. The method used         
c     is a power-law approximation for y(x) between any two points .                
c                                                          [Z.I., Mar. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER i, N, N1, N2                                                      
      DOUBLE PRECISION x(N), y(N), integral, pow, C, delint                     

c     set integral to 0 and accumulate result in the loop                       
      integral = 0.0                                                            
c     calculate weight, wgth, and integrate in the same loop                    
      IF (N2.GT.N1) THEN                                                        
        DO i = N1, N2-1                                                         
          pow = dlog(y(i+1)/y(i)) / dlog(x(i+1)/x(i))                           
          C = y(i) / x(i)**pow                                                  
          delint = (x(i+1)**(pow+1)-x(i)**(pow+1.))*C/(pow+1.)                  
c         add contribution to the integral                                      
          integral = integral + delint                                          
        END DO                                                                  
        ELSE                                                                    
c        integral = 0.0                                                         
        integral = y(1)                                                         
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE DoProduct(NN,Yt,pt,p0,j,Prd)                                     
c     =======================================================================       
c     This is an auxiliary subroutine which evaluates a messy expression            
c     needed to calculate normalization constants for a broken power law            
c     density.                                             [Z.I., Aug. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER NN, i, j                                                          
      DOUBLE PRECISION Yt(NN), pt(NN), Prd, p0                                  

      Prd = Yt(1)**(pt(1) - p0)                                                 
      IF (j.GT.1) THEN                                                          
        DO i = 2, j                                                             
          Prd = Prd * Yt(i)**(pt(i) - pt(i-1))                                  
        END DO                                                                  
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE ROMBERG2(a,b,ss8)                                              
c     =======================================================================       
c     This subroutine performs Romberg integration of 8 functions calculated        
c     in trapzd2 (by calling subroutine TWOFUN) on interval [a,b].                  
c     The results are returned in ss8(1..8). Desired accuracy accRomb is            
c     user supplied and comes through COMMON /numerics/ read in from                
c     'numerics.inc'. This subroutine is based on slightly changed versions         
c     of 'qromb' and 'qromo' from Numerical Recipes.                                
c                                                        [MN & ZI,Aug'96]           
c     =======================================================================       
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER fconv(8),JMAX,JMAXP,K,KM, J, iC, idone, kaux                      
      PARAMETER (JMAX=50, JMAXP=JMAX+1, K=5, KM=K-1)                            
      DOUBLE PRECISION ss, ss8(8), S2D(8,JMAXP), h(JMAXP), sjKM(JMAXP),         
     &                 a, b, EPS, h0, dss, s8(8), chk(8)                        

      EPS = accRomb                                                             
      h0 = 0.0                                                                  
      h(1)=1.0                                                                  
c     intialize convergence flags                                               
      DO iC = 1, 8                                                              
         fconv(iC) = 0                                                          
      END DO                                                                    
c     integrate until all 8 intergrals converge                                 
      idone = 0                                                                 
      j = 0                                                                     
      DO WHILE(idone.NE.1.and.j.LE.JMAX)                                        
        j = j + 1                                                               
c       integrate with j division points                                        
        call trapzd2(a,b,s8,j)                                                  
        DO iC = 1, 8                                                            
           S2D(iC,j) = S8(iC)                                                   
        END DO                                                                  
c       check if any of 8 integrals has converged                               
        IF (j.ge.K) THEN                                                        
           idone = 1                                                            
           DO iC = 1, 8                                                         
             IF (fconv(iC).EQ.0) THEN                                           
c              generate array for polint                                        
               DO kaux = 1, j                                                   
                 sjKM(kaux) = S2D(iC,kaux)                                      
               END DO                                                           
c              predict the integral for stepsize h->h0=0.0                      
               CALL polint(h(j-KM),sjKM(j-KM),K,h0,ss,dss)                      
               IF (dabs(dss).le.EPS*dabs(ss)) THEN                              
                 SS8(iC) = ss                                                   
                 fconv(iC) = 1                                                  
               ELSE                                                             
                 chk(iC) = dabs(dss)/dabs(ss)                                   
               END IF                                                           
             END IF                                                             
             idone = idone*fconv(iC)                                            
           END DO                                                               
        END IF                                                                  
        h(j+1)=0.25*h(j)                                                        
      END DO                                                                    
      IF (j.GE.jMAX) THEN                                                       
        write(*,*)' Reached the limiting number of steps in ROMBERG2'           
        write(*,*)'You might want to change accRomb in the input file'          
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE ROMBY(fnc,a,b,ss)                                              
c     =======================================================================       
c     This subroutine performs Romberg integration of function func on              
c     interval [a,b]. The result is returned in ss. Desired accuracy is set         
c     to 0.002.                                            [Z.I., Feb. 1996]        
c     =======================================================================       
      INTEGER JMAX,JMAXP,K,KM, J                                                
      PARAMETER (JMAX=30, JMAXP=JMAX+1, K=3, KM=K-1)                            
      DOUBLE PRECISION a,b,fnc,ss,EPS, aux, dss,h(JMAXP),s(JMAXP)               
      EXTERNAL fnc                                                              

      EPS = 0.002                                                               
      h(1)=1.                                                                   
      do 11 j=1,JMAX                                                            
        call trapzd(fnc,a,b,s(j),j)                                             
        if (j.ge.K) then                                                        
          aux = 0.0                                                             
          call polint(h(j-KM),s(j-KM),K,aux,ss,dss)                             
          IF (dabs(dss).le.EPS*dabs(ss)) RETURN                                 
        endif                                                                   
        s(j+1)=s(j)                                                             
        h(j+1)=0.25*h(j)                                                        
11    continue                                                                  

      RETURN                                                                    
      END                                                                       


      SUBROUTINE ScaleTo1(Nmax,N,Y)                                             
c     =======================================================================       
c     This subroutine scales vector Y such that Y(1) = 1.0                          
c                                                          [Z.I., Jan. 1997]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER Nmax, N, i                                                        
      DOUBLE PRECISION Y(Nmax), Scale                                           

      Scale = Y(1)                                                              
      DO i = 1, N                                                               
        Y(i) = Y(i) / Scale                                                     
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE SHIFT(X,Nmax,N,Xins,i)                                         
c     =======================================================================       
c     Rearranges a vector X by inserting a new element Xins.    [MN, Aug'96]        
c     =======================================================================       
      implicit none                                                             
      integer Nmax, N, i,j                                                      
      DOUBLE PRECISION X(Nmax),Xins                                             

      DO j = N+1, i+2, -1                                                       
        x(j) = x(j-1)                                                           
      END DO                                                                    
      x(i+1) = xins                                                             

      RETURN                                                                    
      END                                                                       


      SUBROUTINE SIMPSON(N,N1,N2,x,y,integral)                                  
c     =======================================================================       
c     This subroutine calculates integral I(y(x)*dx). Both y and x are              
c     1D arrays, y(i), x(i) with i=1,N (declared with NN). Lower and upper          
c     integration limits are x(N1) and x(N2), respectively. The method used         
c     is Simpson (trapezoid) approximation. The resulting integral is sum of        
c     y(i)*wgth, i=N1,N2.                                  [Z.I., Mar. 1996]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER i, N, N1, N2                                                      
      DOUBLE PRECISION x(N), y(N), wgth, integral                               

c     set integral to 0 and accumulate result in the loop                       
      integral = 0.0                                                            
c     calculate weight, wgth, and integrate in the same loop                    
      IF (N2.GT.N1) THEN                                                        
        DO i = N1, N2                                                           
c         weigths                                                               
          IF (i.NE.N1.AND.i.NE.N2) THEN                                         
            wgth = 0.5 * (x(i+1)-x(i-1))                                        
          ELSE                                                                  
            IF (i.eq.N1) wgth = 0.5 * (x(N1+1)-x(N1))                           
            IF (i.eq.N2) wgth = 0.5 * (x(N2)-x(N2-1))                           
          END IF                                                                
c         add contribution to the integral                                      
          integral = integral + y(i) * wgth                                     
        END DO                                                                  
      ELSE                                                                      
        integral = 0.0                                                          
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE SORT(RA,N)                                                     
c     =======================================================================       
      INTEGER N                                                                 
      DOUBLE PRECISION RA(N)                                                    

      L=N/2+1                                                                   
      IR=N                                                                      
10    CONTINUE                                                                  
        IF(L.GT.1)THEN                                                          
          L=L-1                                                                 
          RRA=RA(L)                                                             
        ELSE                                                                    
          RRA=RA(IR)                                                            
          RA(IR)=RA(1)                                                          
          IR=IR-1                                                               
          IF(IR.EQ.1)THEN                                                       
            RA(1)=RRA                                                           
            RETURN                                                              
          ENDIF                                                                 
        ENDIF                                                                   
        I=L                                                                     
        J=L+L                                                                   
20      IF(J.LE.IR)THEN                                                         
          IF(J.LT.IR)THEN                                                       
            IF(RA(J).LT.RA(J+1))J=J+1                                           
          ENDIF                                                                 
          IF(RRA.LT.RA(J))THEN                                                  
            RA(I)=RA(J)                                                         
            I=J                                                                 
            J=J+J                                                               
          ELSE                                                                  
            J=IR+1                                                              
          ENDIF                                                                 
        GO TO 20                                                                
        ENDIF                                                                   
        RA(I)=RRA                                                               
      GO TO 10                                                                  

      END                                                                       


      SUBROUTINE Spline(x,y,n,yp1,ypn,y2)                                       
c     =======================================================================       
      INTEGER n,NMAX                                                            
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)                                  
      PARAMETER (NMAX=500)                                                      
      INTEGER i,k                                                               
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)                                      

      if (yp1.gt..99e30) then                                                   
        y2(1)=0.                                                                
        u(1)=0.                                                                 
      else                                                                      
        y2(1)=-0.5                                                              
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)                     
      endif                                                                     
      do 11 i=2,n-1                                                             
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))                                       
        p=sig*y2(i-1)+2.                                                        
        y2(i)=(sig-1.)/p                                                        
        u(i)=(6.*((y(i+1)-y(i))/(x(i+                                           
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*                
     *u(i-1))/p                                                                 
11    continue                                                                  
      if (ypn.gt..99e30) then                                                   
        qn=0.                                                                   
        un=0.                                                                   
      else                                                                      
        qn=0.5                                                                  
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))                 
      endif                                                                     
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)                                      
      do 12 k=n-1,1,-1                                                          
        y2(k)=y2(k)*y2(k+1)+u(k)                                                
12    continue                                                                  

      RETURN                                                                    
      END                                                                       


      SUBROUTINE SPLINE2(x,fun,N,coef)                                          
c     =======================================================================       
c     This subroutine finds coefficients coef(i,j) such that                        
c     fun(x)=coef(i,1) + coef(i,2)*x + coef(i,3)*x^2 + coef(i,4)*x^3                
c     for x(i).LE.x.LE.x(i+1) is a cubic spline approximation of fun(x),            
c     with i=1..N.                                         [Z.I., Feb. 1995]        
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER N, i                                                              
      DOUBLE PRECISION x(npY), coef(npY,4), secnder(npY), y2at1, y2atN,         
     &       Dd, xL, xR, dR, dL, fun(npY), fL, fR                               

c     find second derivative, secnder                                           
        y2at1 = (fun(2)-fun(1))/(x(2)-x(1))                                     
        y2atN = (fun(N)-fun(N-1))/(x(N)-x(N-1))                                 
        CALL SPLINE(x,fun,N,y2at1,y2atN,secnder)                                
c     generate coef(i,j), j=1,2,3,4                                             
      DO i = 1, N-1                                                             
        Dd = x(i+1) - x(i)                                                      
        xL = x(i)                                                               
        xR = x(i+1)                                                             
        dL = secnder(i)                                                         
        dR = secnder(i+1)                                                       
        fL = fun(i)                                                             
        fR = fun(i+1)                                                           
        coef(i,1) = (xR*fL-xL*fR)/Dd + dL*xR*Dd/6.*((xR/Dd)**2.-1.)             
        coef(i,1) = coef(i,1) - dR*xL*Dd/6. *((xL/Dd)**2.-1.)                   
        coef(i,2) = (fR-fL)/Dd + dL*Dd/6.*(1.-3.*(xR/Dd)**2.)                   
        coef(i,2) = coef(i,2) - dR*Dd/6.*(1.-3.*(xL/Dd)**2.)                    
        coef(i,3) = (dL*xR-dR*xL)/Dd/2.                                         
        coef(i,4) = (dR-dL)/6./Dd                                               
      END DO                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE splint(xa,ya,y2a,n,x,y)                                        
c     =======================================================================       
      INTEGER n                                                                 
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)                                   
      INTEGER k,khi,klo                                                         
      DOUBLE PRECISION a,b,h                                                    

      klo=1                                                                     
      khi=n                                                                     
1     if (khi-klo.gt.1) then                                                    
        k=(khi+klo)/2                                                           
        if(xa(k).gt.x)then                                                      
          khi=k                                                                 
        else                                                                    
          klo=k                                                                 
        endif                                                                   
      goto 1                                                                    
      endif                                                                     
      h=xa(khi)-xa(klo)                                                         
      if (h.eq.0.) pause 'bad xa input in splint'                               
      a=(xa(khi)-x)/h                                                           
      b=(x-xa(klo))/h                                                           
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**          
     *2)/6.                                                                     

      return                                                                    
      END                                                                       


      SUBROUTINE trapzd(func,a,b,s,n)                                           
c     =======================================================================       
      INTEGER n                                                                 
      DOUBLE PRECISION a,b,s,func                                               
      EXTERNAL func                                                             
      INTEGER it,j                                                              
      DOUBLE PRECISION del,sum,tnm,x                                            

      IF (n.eq.1) THEN                                                          
        s=0.5*(b-a)*(func(a)+func(b))                                           
      ELSE                                                                      
        it=2**(n-2)                                                             
        tnm=it                                                                  
        del=(b-a)/tnm                                                           
        x=a+0.5*del                                                             
        sum=0.                                                                  
        DO j = 1, it                                                            
          sum=sum+func(x)                                                       
          x=x+del                                                               
        END DO                                                                  
        s=0.5*(s+(b-a)*sum/tnm)                                                 
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE WriteOut(is,nG,nameQ,nameNK)                                   
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
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER is, iG, nG, i, length                                             
      CHARACTER*72 strpow, aux, src, chaux*3                                    
      CHARACTER*(*) nameQ(npG), nameNK(10)                                      

      IF (denstyp.eq.0) THEN                                                    
       IF (is.eq.1) THEN                                                        
        src = 'Left-side source spectrum described by'                          
       ELSE                                                                     
        src = 'Right-side source spectrum described by'                         
       END IF                                                                   
      ELSE                                                                      
        src = 'Central source spectrum described by'                            
      END IF                                                                    
      CALL Clean(src, aux, length)                                              
                                                                                
c     #1: black body(ies) for startyp=1                                         
      IF (startyp(is).EQ.1) THEN                                                
       IF (nBB(is).GT.1) THEN                                                   
         CALL ATTACH(aux, length, ' ', src)                                     
c        multiple black bodies                                                  
         write(12,'(a2,a37,i2,a13)')'  ', src, nBB(is),' black bodies'          
         write(12,'(a27)')' with temperatures (in K):'                          
         write(12,'(2x,1p,10e10.3)')(Tbb(is,i),i=1,nBB(is))                     
         write(12,'(a42)')' and relative luminosities, respectively:'           
         write(12,'(1p,10e10.1)')(rellum(is,i),i=1,nBB(is))                     
       ELSE                                                                     
c       for a single black body:                                                
        CALL ATTACH(aux,length,' a black body',src)                             
         write(12,'(a2,a)') '  ',src                                            
         IF (Tstar.LT.9999.999) THEN                                            
           CALL getfs(Tstar,0,1,strpow)                                         
           write(12,'(a19,a5,a2)')' with temperature:',strpow,' K'              
         ELSE                                                                   
           CALL getfs(Tstar,0,1,strpow)                                         
           write(12,'(a19,a6,a2)')' with temperature:',strpow,' K'              
         END IF                                                                 
       END IF                                                                   
      END IF                                                                    
                                                                                
c     #2: Engelke-Marengo function for startyp=2                                
      IF (startyp(is).EQ.2) THEN                                                
         CALL ATTACH(aux, length,' Engelke-Marengo function', src)              
         write(12,'(a2,a)') '  ',src                                            
         CALL getfs(Tbb(is,1),0,1,strpow)                                       
         write(12,'(a13,a5,a16)')' with Teff =',strpow,' K and depth of'        
         write(12,'(a30,F6.1,a2)')' the SiO absorption feature =',              
     &                              xSiO,' %'                                   
      END IF                                                                    
                                                                                
c     #3: power-law(s) for startyp=3                                            
      IF (startyp(is).EQ.3) THEN                                                
       IF (Nlamtr(is).GT.0) THEN                                                
         CALL ATTACH(aux,length,' power law:',src)                              
         write(12,'(a2,a)') '  ',src                                            
         write(12,*)'    lambda      k'                                         
         DO i = 1, Nlamtr(is)                                                   
           write(12,'(1x,1p,e10.3)')lamtr(is,i)                                 
           write(12,'(11x,1p,e10.3)')klam(is,i)                                 
         END DO                                                                 
         write(12,'(1x,1p,e10.3)')lamtr(is,Nlamtr(is)+1)                        
       ELSE                                                                     
         write(12,*)                                                            
     &      ' Input data for the source spectrum is not good.'                  
         write(12,*)' Changed to a 10000 K black body'                          
       END IF                                                                   
      END IF                                                                    
                                                                                
c     spectrum from a file for startyp=4,5,6                                    
      IF (startyp(is).GE.4.AND.startyp(is).LE.6) THEN                           
        write(12,*)' Stellar spectrum supplied from file:'                      
        write(12,'(a2,a70)') '  ',nameStar(is)                                  
      END IF                                                                    
      IF(is.eq.1)                                                               
     &  write(12,*)' --------------------------------------------'              
                                                                                
      IF(is.eq.1) THEN                                                          
c      2) DUST PROPERTIES                                                       
c      2.1 Chemical Composition                                                 
        write(12,*)' Abundances for supported grains:'                          
        write(12,*)' Sil-Ow Sil-Oc Sil-DL grf-DL amC-Hn SiC-Pg'                 
        write(12,'(6f7.3)')(xC(i),i=1,3),xC(4)+xC(5),(xC(i),i=6,7)              
        IF (top.EQ.2) THEN                                                      
          write(12,*)' Abundances for user supplied grains:'                    
          write(12,'(i6,9i7)')(i,i=1,Nfiles)                                    
          write(12,'(10f7.3)')(xCuser(i),i=1,Nfiles)                            
          write(12,*)' User supplied n and k from:'                             
          DO i = 1, Nfiles                                                      
            write(12,'(a2,i1,a2,a70)')'  ',i,') ',nameNK(i)                     
          END DO                                                                
        END IF                                                                  
c      user supplied cross-sections:                                            
       IF (top.EQ.3) THEN                                                       
        DO iG = 1, nG                                                           
          write(12,*)' Optical properties from file:'                           
          write(12,'(a2,a70)')'  ',nameQ(iG)                                    
        END DO                                                                  
       END IF                                                                   
c      2.2 Grain size distribution                                              
       IF (top.NE.3) THEN                                                       
         IF (szds.EQ.3) THEN                                                    
          chaux = 'KMH'                                                         
         ELSE                                                                   
          chaux = 'MRN'                                                         
         END IF                                                                 
         write(12,'(a2,a3,a19)')'  ',chaux,'size distribution:'                 
         CALL getfs(qsd,1,0,strpow)                                             
         write(12,'(a15,a5)')'      Power q:',strpow                            
         write(12,'(a15,1p,e9.2,a8)')                                           
     &                         ' Minimal size:',a1,' microns'                   
         IF (szds.EQ.3) THEN                                                    
            write(12,'(a22,1p,e9.2,a8)')                                        
     &                            ' Characteristic size:',a2,' microns'         
         ELSE                                                                   
           write(12,'(a15,1p,e9.2,a8)')' Maximal size:',a2,' microns'           
         END IF                                                                 
       END IF                                                                   
       write(12,*)' --------------------------------------------'               
c      2.3 Dust temperature on inner boundary                                   
       DO iG = 1, nG                                                            
        CALL getfs(Tsub(iG),0,1,strpow)                                         
        IF (denstyp.eq.0) THEN                                                  
          write(12,'(a45,a5,a2)')                                               
     &      ' Dust temperature on the slab left boundary:', strpow,' K'         
        ELSE                                                                    
          CALL getfs(Tsub(iG),0,1,strpow)                                       
          write(12,'(a41,a5,a2)')                                               
     &      ' Dust temperature on the inner boundary:', strpow,' K'             
        END IF                                                                  
       END DO                                                                   
       write(12,*)' --------------------------------------------'               
      END IF                                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE zbrac(func,x1,x2,Ntry,succes)                                  
c     =======================================================================       
      INTEGER NTRY, succes                                                      
      DOUBLE PRECISION x1,x2,func,FACTOR                                        
      EXTERNAL func                                                             
      PARAMETER (FACTOR=1.6)                                                    
      INTEGER j                                                                 
      DOUBLE PRECISION f1,f2                                                    

      IF(x1.eq.x2) PAUSE 'you have to guess an initial range in zbrac'          
      f1=func(x1)                                                               
      f2=func(x2)                                                               
      succes=1                                                                  
      DO j = 1, NTRY                                                            
        IF(f1*f2.lt.0.) RETURN                                                  
        IF(DABS(f1).lt.DABS(f2)) THEN                                           
          x1=x1+FACTOR*(x1-x2)                                                  
c         f1=func(x1)                                                           
c     IF's added to prevent breaking of slab case on DEC [MN,Jun'99]            
          IF (x1.LT.(-100.)) THEN                                               
            f1 = 5.0e5                                                          
          ELSE                                                                  
            f1=func(x1)                                                         
          END IF                                                                
        ELSE                                                                    
          x2=x2+FACTOR*(x2-x1)                                                  
          f2=func(x2)                                                           
        END IF                                                                  
      END DO                                                                    
      succes=0                                                                  

      RETURN                                                                    
      END                                                                       


      DOUBLE PRECISION FUNCTION zriddr(func,x1,x2,MAXIT,xacc)                   
c     =======================================================================       
      INTEGER MAXIT                                                             
      DOUBLE PRECISION x1,x2,xacc,func,UNUSED                                   
      PARAMETER (UNUSED=-1.11E30)                                               
      EXTERNAL func                                                             
      INTEGER j                                                                 
c     aux will be used as argument for SIGN function few lines below.           
c     It needs to be REAL to conform to FORTRAN90                               
      REAL aux                                                                  
      DOUBLE PRECISION fh,fl,fm,fnew,s,xh,xl,xm,xnew                            

c     IF's added to prevent breaking of slab case on DEC [MN,Jun'99]            
      IF (x1.LT.(-100.)) THEN                                                   
        fl = 5.0e5                                                              
      ELSE                                                                      
        fl=func(x1)                                                             
      END IF                                                                    
      fh=func(x2)                                                               
      if((fl.gt.0..and.fh.lt.0.).or.(fl.lt.0..and.fh.gt.0.))then                
        xl=x1                                                                   
        xh=x2                                                                   
        zriddr=UNUSED                                                           
        do 11 j=1,MAXIT                                                         
          xm=0.5*(xl+xh)                                                        
          fm=func(xm)                                                           
          s=sqrt(fm**2.-fl*fh)                                                  
          if(s.eq.0.)return                                                     
          aux = fl-fh                                                           
          xnew=xm+(xm-xl)*(sign(1.,aux)*fm/s)                                   
          if (DABS(xnew-zriddr).le.xacc) return                                 
          zriddr=xnew                                                           
          fnew=func(zriddr)                                                     
          if (fnew.eq.0.) return                                                
          if(dsign(fm,fnew).ne.fm) then                                         
            xl=xm                                                               
            fl=fm                                                               
            xh=zriddr                                                           
            fh=fnew                                                             
          else if(dsign(fl,fnew).ne.fl) then                                    
            xh=zriddr                                                           
            fh=fnew                                                             
          else if(dsign(fh,fnew).ne.fh) then                                    
            xl=zriddr                                                           
            fl=fnew                                                             
          else                                                                  
            pause 'never get here in zriddr'                                    
          endif                                                                 
          if(dabs(xh-xl).le.xacc) return                                        
11      continue                                                                
        pause 'zriddr exceed maximum iterations'                                
      else if (fl.eq.0.) then                                                   
        zriddr=x1                                                               
      else if (fh.eq.0.) then                                                   
        zriddr=x2                                                               
      else                                                                      
        pause 'root must be bracketed in zriddr'                                
      endif                                                                     

      RETURN                                                                    
      END                                                                       
