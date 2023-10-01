c =======================================================================       
c This is the include file with I/O subroutines              [MN, Apr'99]       
c =======================================================================       
C     Table of Contents                                                         
C                                                                               
C     CLLOSE                                                                    
C     INPSTAR                                                                   
C     INPUT                                                                     
C     OPPEN                                                                     
C     PROUT                                                                     
C     RDINP                                                                     
C     VAL                                                                       
c =======================================================================       
                                                                                

      SUBROUTINE CLLOSE(error,model,Nmodel)                                     
c     =======================================================================       
c     This subroutine closes output files.             [ZI,Feb'96; MN,Apr'99]       
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
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         
      CHARACTER*72 su1, su2, s3, s4                                             
      INTEGER  error, model, Nmodel                                             
c  -----------------------------------------------------------------------       
c     Close the default output file:                                            
      IF (iERROR.NE.0.AND.iOUT.EQ.1) THEN                                       
        write(12,'(a42,i4)')                                                    
     &                 ' There are some error messages for model:',model        
        write(12,*)                                                             
     &         ' Please check m## file (if not produced then rerun)'            
      END IF                                                                    
      IF (iWARNING.NE.0.AND.iOUT.EQ.1.AND.iERROR.EQ.0) THEN                     
        write(12,'(a36,i4)')                                                    
     &  ' There are some warnings for model:',model                             
        write(12,*)                                                             
     &         ' Please check m## file (if not produced then rerun)'            
      END IF                                                                    
      iCUMM = iCUMM + iERROR + iWARNING                                         
      IF (model.EQ.Nmodel.OR.error.EQ.3.OR.error.EQ.4) THEN                     
      IF (error.NE.3) THEN                                                      
       su1=' =========================================================='        
       su2='============================'                                       
        IF(denstyp.eq.4.OR.RDW) THEN                                            
          write(12,'(a59,a28)')su1,su2                                          
        ELSE                                                                    
         IF(denstyp.NE.0) THEN                                                  
           write(12,'(a59)')su1                                                 
         ELSE                                                                   
          su1=' ====================================================='          
          IF(ksi.GT.0) THEN                                                     
           su2='=============='                                                 
           write(12,'(a54,a14)')su1,su2                                         
          ELSE                                                                  
           su2='====='                                                          
           write(12,'(a54,a5)')su1,su2                                          
          END IF                                                                
         END IF                                                                 
        END IF                                                                  
        write(12,'(a23,1p,e8.1,a8)')'   (1) Optical depth at',lamfid,           
     &                              ' microns'                                  
        IF(denstyp.eq.0) THEN                                                   
c  ----------  for slab output ----------------------------                     
         write(12,*)                                                            
     &'  (2) Bol.flux of the left-side source at the slab left boundary'        
         write(12,*)                                                            
     & '  (3) f1=F/Fe1, where F is the overall bol.flux in the slab'            
         write(12,*)                                                            
     & '  (4) Position of the left slab boundary for L=1E4 Lsun'                
         write(12,*)                                                            
     & '  (5) Dust temperature at the right slab face'                          
         write(12,*)                                                            
     & '  (6) Effective temperature of the left source (in K)'                  
         IF(ksi.GT.0) THEN                                                      
           write(12,*)                                                          
     & '  (7) Effective temperature of the right source (in K)'                 
           write(12,*)'  (8) Maximum error in flux conservation (%)'            
         ELSE                                                                   
           write(12,*)'  (7) Maximum error in flux conservation (%)'            
         END IF                                                                 
        ELSE                                                                    
c    ---------- for spherical shell ----------------------------                
         write(12,*)'  (2) Bolometric flux at the inner radius '                
         write(12,*)'  (3) Inner radius for L=1E4 Lsun'                         
         write(12,*)'  (4) Ratio of the inner to the stellar radius'            
         write(12,*)'  (5) Angular size (in arcsec) when Fbol=1E-6 W/m2'        
         write(12,*)'  (6) Dust temperature at the outer edge (in K)'           
         write(12,*)'  (7) Maximum error in flux conservation (%)'              
         IF(denstyp.eq.4.OR.RDW) THEN                                           
          write(12,*)'  (8) Mass-loss rate (in Msun/yr)'                        
          write(12,*)'  (9) Terminal outflow velocity (in km/s)'                
          write(12,*)'  (10) Upper limit of the stellar mass (in Msun)'         
         END IF                                                                 
        END IF                                                                  
        write(12,*)'================================================='          
        IF(iCUMM.EQ.0) write(12,*)' Everything is OK for all models'            
        IF(iSPP.NE.0) THEN                                                      
         IF(denstyp.eq.0.AND.iSPP.eq.3) THEN                                    
          write(12,*)                                                           
     & ' Tables with spectral properties are in files *.spp and *.zpp'          
         ELSE                                                                   
          write(12,*)' Table with spectral properties is in file *.spp'         
         END IF                                                                 
        END IF                                                                  
       IF(iA.NE.0) THEN                                                         
        IF (iA.EQ.1) THEN                                                       
          write(12,*)' All spectra are in file *.stb'                           
        ELSE                                                                    
         IF (denstyp.EQ.0.AND.iA.EQ.3) THEN                                     
           write(12,*)' Spectra are in files *.s## and *.z##'                   
         ELSE                                                                   
           write(12,*)' Spectra are in files *.s##'                             
         END IF                                                                 
        END IF                                                                  
       END IF                                                                   
       IF (denstyp.NE.0.AND.iC.NE.0) THEN                                       
        IF (abs(iC).EQ.1) THEN                                                  
          write(12,*)' All imaging quantities are in file *.itb'                
        ELSE                                                                    
          write(12,*)' Images are in files *.i##'                               
          IF (iV.NE.0.AND.abs(iC).EQ.3)                                         
     &      write(12,*)' Visibility curves are in files *.v##'                  
        END IF                                                                  
        IF (iC.EQ.-3.AND.iPSF.NE.0) THEN                                        
          write(12,*)' Convolved images are in files *.c##'                     
          IF (psftype.LT.3)                                                     
     &    write(12,*)' Point spread functions are in file *.psf'                
        END IF                                                                  
       END IF                                                                   
       IF(iB.NE.0) THEN                                                         
        IF (iB.EQ.1) THEN                                                       
          write(12,*)' All radial profiles are in file *.rtb'                   
        ELSE                                                                    
          write(12,*)' Radial profiles are in files *.r##'                      
        END IF                                                                  
       END IF                                                                   
       IF (iX.EQ.1)                                                             
     &     write(12,*)' All run-time messages are in file *.mtb'                
       IF (iX.GT.1)                                                             
     &     write(12,*)' Run-time messages are in files *.m##'                   
      ELSE                                                                      
        write(12,*)' Ending calculations for this input file'                   
      END IF                                                                    
      END IF                                                                    
                                                                                
      IF (model.EQ.Nmodel.OR.error.EQ.3.OR.error.EQ.4) THEN                     
        write(12,*)                                                             
     &    '========== THE END =============================='                   
      END IF                                                                    
                                                                                
c     Table with spectral properties                                            
      IF (iSPP.NE.0) THEN                                                       
c       In case of slab: add the zpp table after the spp (if desired)           
        IF(denstyp.eq.0.AND.model.eq.Nmodel.AND.iSPP.NE.3) THEN                 
          s3='###   tau0      Psi      fV       fK       f12    C21  '          
          s4=' C31   C43  b8-13 b14-22 B9.8 B11.4  R9.8-18  '                   
          write(19,'(a49)')                                                     
     &              '# ==============================================='         
          write(19,'(a49)')                                                     
     &              '# Properties of Spectra from the slab left side  '         
          write(19,'(a49)')                                                     
     &              '# -----------------------------------------------'         
          write(19,'(a55,a46)')s3,s4                                            
          DO model = 1, Nmodel                                                  
            write(19,'(a100)') zline(model)                                     
          END DO                                                                
         Close(19)                                                              
         Close(24)                                                              
        END IF                                                                  
      END IF                                                                    
      IF (model.EQ.1.OR.error.EQ.3) THEN                                        
        IF (iPSF.EQ.1.AND.psftype.LT.3) Close(23)                               
      END IF                                                                    
c     conditionally close the spectral files                                    
      IF(iA.EQ.1) THEN                                                          
       IF(model.EQ.Nmodel) close(15)                                            
      ELSE                                                                      
       close(15)                                                                
      END IF                                                                    
c     conditionally close the radial files                                      
      IF(iB.EQ.1) THEN                                                          
       IF(model.EQ.Nmodel) close(16)                                            
      ELSE                                                                      
       close(16)                                                                
      END IF                                                                    
c     conditionally close the imaging files                                     
      IF(abs(iC).EQ.1) THEN                                                     
       IF(model.EQ.Nmodel) close(17)                                            
      ELSE                                                                      
       close(17)                                                                
      END IF                                                                    
      IF(iX.EQ.1) THEN                                                          
       IF(model.EQ.Nmodel) close(18)                                            
      ELSE                                                                      
       close(18)                                                                
      END IF                                                                    
                                                                                
      IF (iPSF.NE.0) Close(21)                                                  
      IF (iV.NE.0) Close(22)                                                    

      RETURN                                                                    
      END                                                                       


      SUBROUTINE InpStar(error,is,nameIn)                                       
c     =======================================================================       
c     This subroutine is for reading the stellar input parameters [MN,Mar'99]       
c     =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      CHARACTER strg*40, nameIn*(*)                                             
      INTEGER error, i, is                                                      
      DOUBLE PRECISION RDINP, sum                                               
      LOGICAL Equal, NoEqual                                                    
      Equal = .True.                                                            
      NoEqual = .False.                                                         

      error = 0                                                                 
c     generic temperature for cases other than BB or Engelke-Marengo shape.     
      Tstar = 10000.0                                                           
c     flag for the external spectrum                                            
      startyp(is) = RDINP(Equal,1)                                              
c     if stellar flag is not btw.1 and 6 stop:                                  
      IF(startyp(is).LT.1 .OR. startyp(is).GT.6) THEN                           
       CALL MSG(11)                                                             
       error = 3                                                                
       goto 999                                                                 
      END IF                                                                    
c     #1: black body(ies) for startyp=1                                         
      IF (startyp(is).EQ.1) THEN                                                
c       number of black bodies                                                  
        nBB(is) = RDINP(Equal,1)                                                
c       stellar temperature(s)                                                  
        Tbb(is,1) = RDINP(Equal,1)                                              
        IF (nBB(is).GT.1) THEN                                                  
c       multiple black bodies                                                   
          Tstar = 0.0                                                           
          DO i = 2, nBB(is)                                                     
            Tbb(is,i) = RDINP(NoEqual,1)                                        
            IF (Tbb(is,i).LE.0.0) THEN                                          
              CALL MSG(8)                                                       
              error = 1                                                         
              GOTO 999                                                          
            END IF                                                              
          END DO                                                                
c         read in relative luminosities                                         
          rellum(is,1) = RDINP(Equal,1)                                         
          sum = rellum(is,1)                                                    
          DO i = 2, nBB(is)                                                     
            rellum(is,i) = RDINP(NoEqual,1)                                     
            sum = sum + rellum(is,i)                                            
          END DO                                                                
          IF (sum.LE.0.0) THEN                                                  
             CALL MSG(7)                                                        
             error = 1                                                          
             GOTO 999                                                           
          END IF                                                                
c         normalize to 1                                                        
          DO i = 1, nBB(is)                                                     
            rellum(is,i) = rellum(is,i) / sum                                   
            Tstar = Tstar + rellum(is,i)*Tbb(is,i)**4                           
          END DO                                                                
          Tstar= Tstar**0.25                                                    
        ELSE                                                                    
c         this is for a single black body - Tbb(is,1)                           
          IF(is.eq.1) THEN                                                      
           Tstar = Tbb(1,1)                                                     
           rellum(1,1) = 1.0                                                    
          ELSE                                                                  
           Tstar = Tbb(2,1)                                                     
           rellum(2,1) = 1.0                                                    
          END IF                                                                
          IF (Tbb(is,1).LE.0.0) THEN                                            
             CALL MSG(8)                                                        
             error = 1                                                          
             GOTO 999                                                           
          ELSE                                                                  
c        end if for first or second source                                      
         END IF                                                                 
c       end if for one or many BB                                               
        END IF                                                                  
c     end if for BB-type                                                        
      END IF                                                                    
                                                                                
c     #2: Engelke-Marengo function for startyp=2                                
      IF (startyp(is).EQ.2) THEN                                                
c        effective stellar temperature                                          
         Tbb(is,1) = RDINP(Equal,1)                                             
         Tstar = Tbb(is,1)                                                      
c        depth of SiO abs.feature in %                                          
         xSiO = RDINP(Equal,1)                                                  
         IF (xSiO.LE.0.0) xSiO = 0.0001                                         
         IF (xSiO.GT.100.0) xSiO = 100.0                                        
      END IF                                                                    
                                                                                
c     #3: power-law(s) for startyp=3                                            
      IF (startyp(is).EQ.3) THEN                                                
c       number of transitions                                                   
        Nlamtr(is)= RDINP(Equal,1)                                              
        IF (Nlamtr(is).GT.0) THEN                                               
          lamtr(is,1) = RDINP(Equal,1)                                          
          DO i = 2, Nlamtr(is)+1                                                
            lamtr(is,i) = RDINP(NoEqual,1)                                      
            IF (lamtr(is,i).LT.lamtr(is,i-1)) THEN                              
             CALL MSG(6)                                                        
             error = 1                                                          
             GOTO 999                                                           
            END IF                                                              
          END DO                                                                
          klam(is,1) = RDINP(Equal,1)                                           
          IF (Nlamtr(is).GT.1) THEN                                             
            DO i = 2, Nlamtr(is)                                                
              klam(is,i) = RDINP(NoEqual,1)                                     
            END DO                                                              
          END IF                                                                
        ELSE                                                                    
          startyp(is) = 1                                                       
          Tstar = 10000.0                                                       
        END IF                                                                  
      END IF                                                                    
                                                                                
c     spectrum from a file for startyp=4,5,6                                    
c     startyp.EQ.4 - file gives lambda*L_lambda                                 
c     startyp.EQ.5 - file gives L_lambda                                        
c     startyp.EQ.6 - file gives Lnu = lambda**2*L_lambda                        
      IF (startyp(is).EQ.4.OR.startyp(is).EQ.5.OR.startyp(is).EQ.6) THEN        
        strg = 'spectral shape of external radiation:'                          
        CALL FileMSG(nameStar(is),strg)                                         
      END IF                                                                    
                                                                                
999   RETURN                                                                    
      END                                                                       
                                                     

      SUBROUTINE INPUT(nameIn,nG,nameOut,nameQ,nameNK,TAU1,TAU2,TAUin,          
     &                 Nrec,GridType,Nmodel,error,version)                      
c     =======================================================================       
c     This subroutine reads input data from the file 'filename'. It                 
c     utilizes the subroutine RDINP written by Moshe Elitzur.                       
c                                                          [Z.I., Nov. 1995]        
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
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         
      CHARACTER LamStr(20)*72, strpow*72, strg*40, version*3                    
      CHARACTER*(*) nameIn, nameOut, nameQ(npG), nameNK(10), namePSF*70,        
     &         nameTAU*70                                                       
      INTEGER iG, nG, Nmodel, i, EtaOK, error, GridType, istop,                 
     &        ioverflw, Nrec, Nmax                                              
c     Nmax is the size of user supplied ETA file                                
      PARAMETER (Nmax = 1000)                                                   
      DOUBLE PRECISION RDINP, TAU1, TAU2, sum, a, b, x(Nmax), e(Nmax),          
     &       aa(Nmax), bb(Nmax), TAUin(Nrec), Ceta, x1, psf1                    
      LOGICAL Equal, NoEqual                                                    
      Equal = .True.                                                            
      NoEqual = .False.                                                         

      error = 0                                                                 
c     open output file                                                          
      open(12,file=nameOut,STATUS='UNKNOWN')                                    
      write(12,*)'==========================='                                  
      write(12,*)' Output from program Dusty '                                  
      write(12,*)' Version: ',version                                           
      write(12,*)'==========================='                                  
      write(12,*)' '                                                            
      write(12,*)' INPUT parameters from file:'                                 
      write(12,'(2x,a70)')nameIn                                                
      write(12,*)' '                                                            
c     open input file                                                           
      open(1,ERR=998,file=nameIn,STATUS='OLD')                                  
      rewind(1)                                                                 
                                                                                
c     *************************                                                 
c     ** PHYSICAL PARAMETERS **                                                 
c     *************************                                                 
                                                                                
c     1) EXTERNAL RADIATION                                                     
      CALL InpStar(error,1,nameIn)                                              
      IF(error.NE.0) goto 996                                                   
                                                                                
c     2) DUST PROPERTIES                                                        
c     # of different dust grains, to be used in a future version                
c     nG = RDINP(Equal,1)                                                       
      nG = 1                                                                    
c     2.1 Chemical Composition                                                  
c     type of optical properties                                                
      top = RDINP(Equal,1)                                                      
      IF (top.NE.1.AND.top.NE.2.AND.top.NE.3) THEN                              
        CALL MSG(9)                                                             
        error = 1                                                               
        GOTO 999                                                                
      END IF                                                                    
c     for top.LT.3 read in abundances for supported grains                      
      IF (top.LT.3) THEN                                                        
        xC(1) = RDINP(Equal,1)                                                  
        IF (xC(1).LT.0.0) xC(1) = 0.0                                           
        sum = xC(1)                                                             
        DO i = 2, 7                                                             
c         special care to be taken of graphite (1/3-2/3 rule):                  
          IF (i.NE.5) THEN                                                      
            xC(i) = RDINP(NoEqual,1)                                            
            IF (xC(i).LT.0.0) xC(i) = 0.0                                       
c           graphite (perpendicular to c axis) :                                
            IF (i.EQ.4) xC(i) = 2.*xC(i)/3.                                     
          ELSE                                                                  
c           graphite (parallel to c axis) :                                     
            xC(i) = 0.5 * xC(i-1)                                               
          END IF                                                                
          sum = sum + xC(i)                                                     
        END DO                                                                  
      END IF                                                                    
c     user supplied n and k:                                                    
      IF (top.EQ.2) THEN                                                        
        Nfiles = RDINP(Equal,1)                                                 
c       file names                                                              
        strg = 'optical constants:'                                             
        DO i = 1, Nfiles                                                        
          CALL FileMSG(nameNK(i),strg)                                          
        END DO                                                                  
        IF(error.NE.0) goto 996                                                 
c       abundances                                                              
        xCuser(1) = RDINP(Equal,1)                                              
        IF (xCuser(1).LT.0.0) xCuser(1) = 0.0                                   
        sum = sum + xCuser(1)                                                   
        IF (Nfiles.GT.1) THEN                                                   
          DO i = 2, Nfiles                                                      
            xCuser(i) = RDINP(NoEqual,1)                                        
            IF (xCuser(i).LT.0.0) xCuser(i) = 0.0                               
            sum = sum + xCuser(i)                                               
          END DO                                                                
        END IF                                                                  
      END IF                                                                    
      IF (top.LT.3) THEN                                                        
        IF (sum.LE.0.0) THEN                                                    
          CALL MSG(5)                                                           
          error = 1                                                             
          GOTO 999                                                              
        END IF                                                                  
c       normalize abundances for supported grains:                              
        DO i = 1, 7                                                             
          xC(i) = xC(i) / sum                                                   
        END DO                                                                  
c       normalize abundances for user supplied grains                           
        IF (top.EQ.2) THEN                                                      
          DO i = 1, Nfiles                                                      
            xCuser(i) = xCuser(i) / sum                                         
          END DO                                                                
        END IF                                                                  
      END IF                                                                    
c     user supplied cross-sections:                                             
      IF (top.EQ.3) THEN                                                        
c       filename for Qabs and Qsca                                              
        strg= 'abs. and scatt. cross-sections:'                                 
        DO iG = 1, nG                                                           
          CALL FileMSG(nameQ(iG),strg)                                          
        END DO                                                                  
      END IF                                                                    
                                                                                
c     2.2 Grain size distribution                                               
      IF (top.NE.3) THEN                                                        
c       type of size distribution                                               
        szds = RDINP(Equal,1)                                                   
        IF (szds.NE.1.AND.szds.NE.2.AND.szds.NE.3) THEN                         
          CALL MSG(10)                                                          
          error = 1                                                             
          GOTO 999                                                              
        END IF                                                                  
c       grain sizes                                                             
        IF (szds.GT.1) THEN                                                     
          qsd = RDINP(Equal,1)                                                  
          a1 = RDINP(Equal,1)                                                   
          IF (a1.LE.0.0) a1 = 0.0001                                            
          a2 = RDINP(Equal,1)                                                   
          IF (szds.EQ.2.AND.a2.LT.a1) a2 = a1                                   
        ELSE                                                                    
          qsd = 3.5                                                             
          a1 = 0.005                                                            
          a2 = 0.25                                                             
        END IF                                                                  
      END IF                                                                    
                                                                                
c     2.3 Dust temperature on inner boundary                                    
      DO iG = 1, nG                                                             
        Tsub(iG) = RDINP(Equal,1)                                               
      END DO                                                                    
                                                                                
c     3) DENSITY DISTRIBUTION                                                   
                                                                                
c     parameter describing eta function:                                        
      denstyp = RDINP(Equal,1)                                                  
c     WriteOut prints all input data, read so far, in fname.out 	               
c     It has to be placed here b/c geometry (spherical or slab) is not determine
c     until denstyp is read [MN,Sep'99].                                        
      CALL WriteOut(1,nG,nameQ,nameNK)                                          
c     The following transfers are caused by the change of order                 
c     of density types [Dec'96]                                                 
      IF (denstyp.EQ.1) THEN                                                    
        denstyp = 2                                                             
        ELSE                                                                    
        IF (denstyp.EQ.2) THEN                                                  
          denstyp = 3                                                           
          ELSE                                                                  
          IF (denstyp.EQ.3) THEN                                                
            denstyp = 5                                                         
            ELSE                                                                
            IF (denstyp.EQ.4) THEN                                              
              denstyp = 4                                                       
              ELSE                                                              
              IF (denstyp.EQ.5) denstyp = 7                                     
            END IF                                                              
          END IF                                                                
        END IF                                                                  
      END IF                                                                    
c     initialize the logical variable RDW (in 'denstyp.inc', used in ChkSplin)  
      IF(denstyp.eq.5.OR.denstyp.eq.6) THEN                                     
        RDW = .true.                                                            
      ELSE                                                                      
        RDW = .false.		  	                                                      
      END IF                                                                    
                                                                                
      EtaOK = 0                                                                 
      Ntr = 0                                                                   
c     Read parameters for each type of density distribution                     
c     slab geometry                                                             
      IF (denstyp.EQ.0) THEN                                                    
        EtaOK = 1                                                               
        write(12,'(a33)') ' Calculation in planar geometry:'                    
        mu1 = RDINP(Equal,1)                                                    
        IF (mu1.GT.1.0) mu1 = 1.0                                               
        IF (mu1.LT.0.0) mu1 = -1.0                                              
        write(12,'(a35,1p,e11.3)')                                              
     &                ' cos of left illumination angle = ',mu1                  
        ksi = RDINP(Equal,1)                                                    
        IF (ksi.LT.0.0) ksi = 0.0                                               
        IF (ksi.GT.1.0) THEN                                                    
          ksi = 1.0                                                             
        END IF                                                                  
        write(12,'(a6,1p,e11.3)')'  R = ',ksi                                   
        IF(ksi.GT.0.) THEN                                                      
c         in case of additional source on the right                             
          mu2 = RDINP(Equal,1)                                                  
          IF (mu2.GT.1.0) mu2 = 1.0                                             
          IF (mu2.LT.0.0) mu2 = -1.0                                            
          write(12,'(a36,1p,e10.3)')                                            
     &                ' cos of right illumination angle = ',mu2                 
          CALL InpStar(error,2,nameIn)                                          
          CALL WriteOut(2,nG,nameQ,nameNK)                                      
        ELSE                                                                    
c         even if no second source is supplied, mu2 needs a value               
c         of 1 to avoid crashing in the formulae in SLBStar [MN]                
          mu2 = 1.                                                              
        END IF                                                                  
      ELSE                                                                      
        ksi = 0.                                                                
      END IF                                                                    
c     smooth or broken power laws                                               
      IF (denstyp.EQ.1.OR.denstyp.EQ.2) THEN                                    
        EtaOK = 1                                                               
        Ntr = RDINP(Equal,1)                                                    
c       changed definition                                                      
        Ntr = Ntr - 1                                                           
c       read in transition radii                                                
        IF (Ntr.GT.0) THEN                                                      
          Ytr(1) = RDINP(Equal,1)                                               
          IF (Ntr.GT.1) THEN                                                    
            DO i = 2, Ntr                                                       
              Ytr(i) = RDINP(NoEqual,1)                                         
            END DO                                                              
          END IF                                                                
          Yout = RDINP(NoEqual,1)                                               
        ELSE                                                                    
          Yout = RDINP(Equal,1)                                                 
        END IF                                                                  
        IF (Yout.LE.1.0) Yout = 1.001                                           
c       read in powers                                                          
        pow = RDINP(Equal,1)                                                    
        IF (Ntr.GT.0) THEN                                                      
          DO i = 1, Ntr                                                         
            ptr(i) = RDINP(NoEqual,1)                                           
          END DO                                                                
        END IF                                                                  
c       print info to the output file                                           
        IF (Ntr.EQ.0) THEN                                                      
          CALL getfs(pow,2,0,strpow)                                            
          write(12,'(a38,a5)')                                                  
     &      ' Density described by 1/r**k with k =',strpow                      
          write(12,'(a21,1p,e10.3)')'  Relative thickness:',Yout                
        ELSE                                                                    
          IF (denstyp.EQ.1) THEN                                                
            write(12,*)' Density described by smooth power law'                 
          ELSE                                                                  
            write(12,*)' Density described by a broken power law:'              
          END IF                                                                
          write(12,*)'  power   Ytransition'                                    
          write(12,*)'  -------------------'                                    
          write(12,*)'              1.0'                                        
          CALL getfs(pow,2,0,strpow)                                            
          write(12,'(a2,a5)')'  ',strpow                                        
          DO i = 1, Ntr                                                         
            write(12,'(a10,1p,e10.3)')'          ',Ytr(i)                       
            CALL getfs(ptr(i),2,0,strpow)                                       
            write(12,'(a2,a5)')'  ',strpow                                      
          END DO                                                                
          write(12,'(a10,1p,e10.3)')'          ',Yout                           
        END IF                                                                  
      END IF                                                                    
c     exponential law                                                           
      IF (denstyp.EQ.3) THEN                                                    
        EtaOK = 1                                                               
        Yout = RDINP(Equal,1)                                                   
        IF (Yout.LE.1.0) Yout = 1.001                                           
        pow = RDINP(Equal,1)                                                    
        IF (pow.LE.0.0) THEN                                                    
          EtaOK = 0                                                             
          ELSE                                                                  
          write(12,*)' Density described by exponential distribution'           
          write(12,'(a21,1p,e10.3)')'               Sigma:',pow                 
          write(12,'(a21,1p,e10.3)')'  Relative thickness:',Yout                
        END IF                                                                  
      END IF                                                                    
c     default approximation and default numerics for rad. driven winds          
      IF (denstyp.EQ.4.OR.denstyp.EQ.5) THEN                                    
        EtaOK = 1                                                               
        Yout = RDINP(Equal,1)                                                   
        IF (Yout.LE.1.0) Yout = 1.001                                           
c       ** DEFAULT ** for epsilon = v1/ve = u1/ue:                              
        pow = 0.2                                                               
        IF(denstyp.EQ.5) THEN                                                   
c         ** DEFAULT ** for Max(GravCor = Fgrav/Frad_press):                    
          ptr(1) = 0.5                                                          
c         convergence criterion:                                                
          ptr(2) = 1.0                                                          
        END IF                                                                  
        write(12,*)' Density for radiatively driven winds from'                 
        IF (denstyp.EQ.4) THEN                                                  
          write(12,*)' analytic approximation for gray dust.'                   
        ELSE                                                                    
          write(12,*)' full dynamic calculation.'                               
        END IF                                                                  
        write(12,'(a21,1p,e10.3)')'  Relative thickness:',Yout                  
      END IF                                                                    
c     full dynamical calculation for radiatively driven winds                   
      IF (denstyp.EQ.6) THEN                                                    
        EtaOK = 1                                                               
c       Yout:                                                                   
        Yout = RDINP(Equal,1)                                                   
        IF (Yout.LE.1.0) Yout = 1.001                                           
c       u1:                                                                     
        pow = RDINP(Equal,1)                                                    
c       GravCor:                                                                
        ptr(1) = RDINP(Equal,1)                                                 
c       convergence criterion:                                                  
        ptr(2) = 1.0                                                            
c       default limits                                                          
        IF (ptr(1).GT.1.0) ptr(1) = 0.9999                                      
        IF (ptr(1).LT.0.0) ptr(1) = 0.0                                         
        write(12,*)' Density for radiatively driven winds:'                     
        write(12,*)' Full dynamic calculation (denstyp.EQ.6),'                  
        write(12,'(a6,1p,e10.3,a20,e10.3)')                                     
     &                   '  u1 =',pow,' Gamma^{-1}_{max} =',ptr(1)              
      END IF                                                                    
c     user specified table for ETA                                              
      IF(denstyp.EQ.7) THEN                                                     
        EtaOK = 1                                                               
        strg = 'dust density distribution:'                                     
        CALL FileMSG(nameETA,strg)                                              
        write(12,*)' Density distribution supplied from file:'                  
        write(12,'(2x,a70)') nameETA                                            
c       read in the density                                                     
        open(31,ERR=997,file=nameETA,STATUS='OLD')                              
c       three lines in the header:                                              
        DO i = 1, 3                                                             
          read(31,*,ERR=997)                                                    
        END DO                                                                  
        istop = 0                                                               
        i = 0                                                                   
        DO WHILE (istop.GE.0)                                                   
          read(31,*,END=900,ERR=997,iostat=istop)a, b                           
          IF (istop.GE.0) THEN                                                  
            i = i + 1                                                           
            x(i) = a                                                            
            e(i) = b                                                            
            IF (i.EQ.1) x1 = x(i)                                               
            yEta7(i) = x(i) / x1                                                
          END IF                                                                
        END DO                                                                  
900     close(31)                                                               
        Eta7OK = 7                                                              
        nYEta7 = i                                                              
        IF (nYEta7.LT.2) goto 997                                               
c       if input positions in descending order turn them around                 
        IF (yEta7(1).GT.yEta7(2)) THEN                                          
          DO i = 1, nYEta7                                                      
            aa(i) = yEta7(i)                                                    
            bb(i) = e(i)                                                        
          END DO                                                                
          DO i = 1, nYEta7                                                      
            yEta7(i) = aa(nYEta7+1-i)                                           
            e(i) = bb(nYEta7+1-i)                                               
          END DO                                                                
        END IF                                                                  
c       relative thickness                                                      
        Yout = yEta7(nYEta7)                                                    
        write(12,'(a21,1p,e10.3)')'  Relative thickness:',Yout                  
        IF (Yout.LE.1.0) Yout = 1.001                                           
c       integrate and ...                                                       
        CALL SIMPSON(Nmax,1,nYEta7,yEta7,e,Ceta)                                
c       ... renormalize                                                         
        DO i = 1, nYEta7                                                        
          Eta7(i) = e(i) / Ceta                                                 
        END DO                                                                  
      END IF                                                                    
c     done with the reading of density distribution                             
      IF (EtaOK.NE.1) THEN                                                      
        CALL MSG(3)                                                             
        error = 1                                                               
        GOTO 999                                                                
      END IF                                                                    
      write(12,*)' --------------------------------------------'                
                                                                                
c     4) OPTICAL DEPTH                                                          
c     grid type                                                                 
      GridType = RDINP(Equal,1)                                                 
      IF (GridType.EQ.3) THEN                                                   
c     TAU-grid from a file                                                      
        strg = 'user supplied TAU-grid:'                                        
        CALL FileMSG(nameTAU,strg)                                              
c       read optical depths                                                     
        open(32,ERR=992,file=nameTAU,STATUS='OLD')                              
c       fiducial wavelength                                                     
c       (the second argument of RDINP is the unit)                              
        lamfid = RDINP(Equal,32)                                                
        i = 0                                                                   
        DO WHILE (i.LT.Nrec)                                                    
          read(32,*,END=902,ERR=992) a                                          
          i = i + 1                                                             
          TAUin(i) = a                                                          
        END DO                                                                  
902     close(32)                                                               
        Nmodel = i                                                              
c       sort the tau-grid                                                       
        CALL Sort(TAUin,Nmodel)                                                 
        TAU1 = TAUin(1)                                                         
        IF (TAU1.LE.0.0) TAU1 = 0.0001                                          
        TAU2 = TAUin(Nmodel)                                                    
      ELSE                                                                      
c      fiducial wavelength                                                      
       lamfid = RDINP(Equal,1)                                                  
c      total optical depths at lamfid                                           
       TAU1 = RDINP(Equal,1)                                                    
       IF (TAU1.LE.0.0) TAU1 = 0.0001                                           
       TAU2 = RDINP(Equal,1)                                                    
       IF (TAU2.LE.TAU1) THEN                                                   
         TAU2 = TAU1                                                            
         Nmodel = 1                                                             
       END IF                                                                   
c      read number of models                                                    
       Nmodel = RDINP(Equal,1)                                                  
       IF (Nmodel.GT.(Nrec-1)) Nmodel = Nrec-1                                  
       IF (Nmodel.LT.1) Nmodel = 1                                              
      END IF                                                                    
      IF (Nmodel.GT.1) THEN                                                     
        write(12,'(a19,1p,e8.1,a8)')' Optical depths at',lamfid,                
     &                              'microns'                                   
        write(12,'(a14,1p,e9.2,a3,e9.2)')' ranging from',TAU1,' to',TAU2        
        IF (GridType.EQ.1) strg=' models with linear grid    '                  
        IF (GridType.EQ.2) strg=' models with logarithmic grid'                 
        IF (GridType.EQ.3) strg=' models with grid from file  '                 
        write(12,'(a1,i4,a)')' ', Nmodel, strg                                  
        IF (GridType.EQ.3) write(12,'(a4,a70)')'    ',nameTAU                   
      ELSE                                                                      
        write(12,'(a18,1p,e8.1,a9,e9.2)')' Optical depth at',lamfid,            
     &                                  ' microns:',TAU1                        
      END IF                                                                    
                                                                                
c     ************************                                                  
c     ** NUMERICAL ACCURACY **                                                  
c     ************************                                                  
c     some of the parameters were originally left to user to                    
c     specify them, most of them are now unchangeable defaults                  
c     max. increase of scaled TAU/TAUtot                                        
c      delTAUsc = RDINP(Equal,1)                                                
c      IF (delTAUsc.LE.0.0) delTAUsc = 0.3                                      
      delTAUsc = 0.3                                                            
c     max. increase in the ratio of two y                                       
c      facc = RDINP(Equal,1)                                                    
c      IF (facc.LE.0.0) facc = 2.0                                              
      facc = 2.0                                                                
c     max allowed ratio of Eta for two consecutive pts.                         
      EtaRat = 4.0                                                              
c     # of rays per radial step                                                 
c      Nins = RDINP(Equal,1)                                                    
c      IF (Nins.LE.0) Nins = 2                                                  
      Nins = 2                                                                  
c     accuracy for numerical integration in ROMBERG                             
c      accRomb = RDINP(Equal,1)                                                 
c      IF (accRomb.LE.0.0) accRomb = 0.0001                                     
c      IF (accRomb.GT.0.005) accRomb = 0.005                                    
      accRomb = 0.0001                                                          
c     accuracy for flux conservation                                            
      accuracy = RDINP(Equal,1)                                                 
      IF (accuracy.LE.0.0) accuracy = 0.05                                      
c     protect against a very large value for accuracies                         
      IF (accuracy.GT.0.25) accuracy = 0.25                                     
c     accuracy for convergence (typical 0.0001)                                 
      accConv = accuracy / 500.0                                                
      accFbol = 10.*accConv                                                     
c     dynamical range                                                           
c      dynrange = RDINP(Equal,1)                                                
c      IF (dynrange.LE.0.0) dynrange = 1.0e-14                                  
c      IF (dynrange.GT.0.0001) dynrange = 1.0e-4                                
      dynrange = 1.0e-15                                                        
      IF (accuracy.GE.0.1) THEN                                                 
        CALL getfs(accuracy*100,0,1,strpow)                                     
        write(12,'(a20,a3,a1)')' Required accuracy:',strpow,'%'                 
      ELSE                                                                      
        CALL getfs(accuracy*100,0,1,strpow)                                     
        write(12,'(a20,a2,a1)')' Required accuracy:',strpow,'%'                 
      END IF                                                                    
      write(12,*)' --------------------------------------------'                
c     ******************                                                        
c     ** OUTPUT FLAGS **                                                        
c     ******************                                                        
c     A flag for additional output for us only [MN]:                            
c     If iInn=1: print err.vs.iter in unt=38 (fname.err) for all models         
c     and additionally list scaled fbol(y) and Ubol(y) in m-files.              
c     iInn replaces iLTS in 'output.inc'                                        
      iInn = 0                                                                  
c     these two by default                                                      
      iSUM = 1                                                                  
      iOUT = 1                                                                  
      iVerb = RDINP(Equal,1)                                                    
c     spectral properties                                                       
      iSPP = RDINP(Equal,1)                                                     
c     spectra                                                                   
      iA = RDINP(Equal,1)                                                       
c     images (intensity) for spherical only                                     
      IF(denstyp.NE.0) THEN                                                     
       iC = RDINP(Equal,1)                                                      
       IF (iC.NE.0) THEN                                                        
        NlambdaOut = RDINP(Equal,1)                                             
        IF (NlambdaOut.GE.1) THEN                                               
          DO i = 1, NlambdaOut                                                  
            LambdaOut(i) = RDINP(NoEqual,1)                                     
c           make sure the wavelengths are inside Dusty's range                  
            IF (LambdaOut(i).LE.0.01) LambdaOut(i) = 0.01                       
            IF (LambdaOut(i).GT.36000) LambdaOut(i) = 36000                     
          END DO                                                                
          ioverflw = 0                                                          
          DO i = 1, NlambdaOut                                                  
            IF (LambdaOut(i).LT.0.995) THEN                                     
              CALL getfs(LambdaOut(i),2,1,LamStr(i))                            
            ELSE                                                                
              IF (LambdaOut(i).LT.9.95) THEN                                    
               CALL getfs(LambdaOut(i),1,0,LamStr(i))                           
              ELSE                                                              
                IF (LambdaOut(i).LT.99.5) THEN                                  
                 CALL getfs(LambdaOut(i),0,0,LamStr(i))                         
                ELSE                                                            
                 CALL getfs(LambdaOut(i),0,1,LamStr(i))                         
c                  IDINT truncates the numbers                                  
c                  prn = IDINT(LambdaOut(i))                                    
c                  write (LamStr(i),'(I6)') prn                                 
                END IF                                                          
              END IF                                                            
            END IF                                                              
            IF (LambdaOut(i).GT.9999.5) THEN                                    
              ioverflw = 1                                                      
              strpow = LamStr(i)                                                
              strpow(4:4) = '*'                                                 
              strpow(5:5) = ' '                                                 
              LamStr(i) = strpow                                                
            END IF                                                              
          END DO                                                                
        END IF                                                                  
        write(12,*)' Images requested for these wavelengths (mic)'              
        write(12,'(a1,20a5)')' ',(LamStr(i),i=1,NlambdaOut)                     
        IF (ioverflw.EQ.1) write(12,*)'  *: in mm'                              
c       convolved images  (only for our use)                                    
        IF (iC.LT.0) THEN                                                       
         iPSF = RDINP(Equal,1)                                                  
         IF (iPSF.NE.0) THEN                                                    
           Theta1 = RDINP(Equal,1)                                              
           write(12,'(a39,1p,e7.1)')                                            
     &           ' Convolved images produced for theta1=',Theta1                
           psftype = RDINP(Equal,1)                                             
           IF (psftype.NE.1.AND.psftype.NE.2.AND.psftype.NE.3) goto 994         
           IF (psftype.LT.3) THEN                                               
c            Gaussians, read in parameters                                      
c            FWHM for the first component                                       
             FWHM1(1) = RDINP(Equal,1)                                          
             IF (NlambdaOut.GT.1) THEN                                          
               DO i = 2, NlambdaOut                                             
                FWHM1(i) = RDINP(NoEqual,1)                                     
               END DO                                                           
             END IF                                                             
             IF (psftype.EQ.2) THEN                                             
c             relative strength for the second component                        
              kPSF(1) = RDINP(Equal,1)                                          
              IF (NlambdaOut.GT.1) THEN                                         
                DO i = 2, NlambdaOut                                            
                  kPSF(i) = RDINP(NoEqual,1)                                    
                END DO                                                          
              END IF                                                            
c             FWHM for the second component                                     
              FWHM2(1) = RDINP(Equal,1)                                         
              IF (NlambdaOut.GT.1) THEN                                         
                DO i = 2, NlambdaOut                                            
                  FWHM2(i) = RDINP(NoEqual,1)                                   
                END DO                                                          
              END IF                                                            
             END IF                                                             
             write(12,*)' The point spread functions are Gaussians'             
           ELSE                                                                 
c           user supplied PSF                                                   
            strg = 'point spread function:'                                     
            CALL FileMSG(namePSF,strg)                                          
            write(12,*)' The point spread function supplied from file'          
            write(12,'(2x,a70)')namePSF                                         
            open(31,ERR=995,file=namePSF,STATUS='OLD')                          
c           three lines in the header:                                          
            DO i = 1, 3                                                         
              read(31,*,ERR=995)                                                
            END DO                                                              
            istop = 0                                                           
            i = 0                                                               
            DO WHILE (istop.GE.0)                                               
              read(31,*,END=901,ERR=995,iostat=istop)a, b                       
              IF (istop.GE.0) THEN                                              
                i = i + 1                                                       
                IF (i.EQ.1) THEN                                                
                  psf1 = b                                                      
                  IF (a.NE.0.0) goto 995                                        
                END IF                                                          
                xpsf(i) = a                                                     
                ypsf(i) = b / psf1                                              
              END IF                                                            
            END DO                                                              
901         close(31)                                                           
            Npsf = i                                                            
c           scale to 1 at the center                                            
            CALL ScaleTo1(1000,Npsf,ypsf)                                       
c           find equivalent FWHM                                                
            istop = 0                                                           
            i = 1                                                               
            DO WHILE (istop.EQ.0)                                               
              i = i + 1                                                         
              IF (ypsf(i).LE.0.5) istop = 1                                     
            END DO                                                              
c           linear interpolation                                                
            FWHM1(1) = (xpsf(i)-xpsf(i-1))/(ypsf(i)-ypsf(i-1))                  
            FWHM1(1) = (FWHM1(1)*(0.5-ypsf(i-1))+xpsf(i-1))*2                   
            FWHM2(1) = 0.0                                                      
            write(12,'(a18,1p,e8.1)')' Equivalent FWHM:',FWHM1(1)               
           END IF                                                               
         END IF                                                                 
        END IF                                                                  
c       visibility (only if the intensity is requested)                         
        iV = RDINP(Equal,1)                                                     
        IF(iV.NE.0) iV = iC                                                     
        write(12,*)' --------------------------------------------'              
       ELSE                                                                     
        iPSF = 0                                                                
        iV = 0                                                                  
       END IF                                                                   
      ELSE                                                                      
       iC = 0                                                                   
      END IF                                                                    
      write(12,*)' '                                                            
c     radial quantities                                                         
      iB = RDINP(Equal,1)                                                       
c     run-time messages                                                         
      iX = RDINP(Equal,1)                                                       
      iINP = 0                                                                  
c     *** DONE ***                                                              
c     if everything is OK, close the input file and finish                      
999   goto 996                                                                  
c     or in the case of err reading files...                                    
992   write(12,*)' ***  FATAL ERROR IN DUSTY  *************'                    
      write(12,*)' File with user supplied TAU-grid:'                           
      write(12,'(2x,a70)') nameTAU                                              
      write(12,*)' is not properly formatted?!'                                 
      write(12,*)' ****************************************'                    
      close(12)                                                                 
      error = 3                                                                 
      goto 996                                                                  
994   CALL MSG(12)                                                              
      close(12)                                                                 
      error = 3                                                                 
      goto 996                                                                  
995   write(12,*)' ***  FATAL ERROR IN DUSTY  *************'                    
      write(12,*)' File with the point spread function:'                        
      write(12,'(a2,a70)')'  ', namePSF                                         
      write(12,*)' is not properly formatted?!'                                 
      write(12,*)' ****************************************'                    
      close(12)                                                                 
      error = 3                                                                 
997   write(12,*)' ***  FATAL ERROR IN DUSTY  *************'                    
      write(12,*)' File with the dust density distribution:'                    
      write(12,'(2x,a70)') nameETA                                              
      write(12,*)' is not properly formatted?!'                                 
      write(12,*)' ****************************************'                    
      close(12)                                                                 
      error = 3                                                                 
998   write(12,*)' ***  FATAL ERROR IN DUSTY  ****'                             
      write(12,*)' Input file:'                                                 
      write(12,'(2x,a70)') nameIn                                               
      write(12,*)' is missing?!'                                                
      write(12,*)' *******************************'                             
      close(12)                                                                 
      error = 3                                                                 

996   close(1)                                                                  
      RETURN                                                                    
      END                                                                       

                                                                                
      SUBROUTINE OPPEN(model,RootName,length)                                   
c     =======================================================================       
c     This subroutine prints the results out.              [Z.I., Feb. 1996]        
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
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         
      CHARACTER ch5*5, RootName*(*), fname*235                                  
      CHARACTER*72 header1, s3, s4                                              
      INTEGER model, length, i                                                  

c     set up the status indicators                                              
      iERROR = 0                                                                
      iWARNING = 0                                                              
      IF (model.EQ.1) iCUMM = 0                                                 
c     The following files pertain to ALL models and are open if model.EQ.1      
      IF (model.EQ.1) THEN                                                      
c      the header to output file *.OUT is moved to PrOut [MN,Sep'99]            
                                                                                
c      open file with spectral properties RootName.SPP                          
       IF (iSPP.NE.0) THEN                                                      
        CALL Attach(RootName,length,'.spp',fname)                               
        open(19,file=fname,STATUS='UNKNOWN')                                    
        IF(denstyp.EQ.0) THEN                                                   
          write(19,'(a49)')                                                     
     &              '# ==============================================='         
          write(19,'(a49)')                                                     
     &              '# Properties of Spectra from the slab right side '         
          write(19,'(a49)')                                                     
     &              '# -----------------------------------------------'         
        ELSE                                                                    
          call line(1,2,19)                                                     
          write(19,'(a23)')'#  Spectral Properties '                            
          call line(1,1,19)                                                     
        END IF                                                                  
         s3='###   tau0      Psi      fV       fK       f12    C21  '           
         s4=' C31   C43  b8-13 b14-22 B9.8 B11.4  R9.8-18  '                    
         write(19,'(a55,a46)')s3,s4                                             
        IF (denstyp.EQ.0.AND.iSPP.EQ.3) THEN                                    
        CALL Attach(RootName,length,'.zpp',fname)                               
          open(24,file=fname,STATUS='UNKNOWN')                                  
          write(24,'(a49)')                                                     
     &              '# ==============================================='         
          write(24,'(a49)')                                                     
     &              '# Properties of Spectra from the slab left side  '         
          write(24,'(a49)')                                                     
     &              '# -----------------------------------------------'         
         s3='###   tau0      Psi      fV       fK       f12    C21  '           
         s4=' C31   C43  b8-13 b14-22 B9.8 B11.4  R9.8-18  '                    
           write(24,'(a55,a46)')s3,s4                                           
        END IF                                                                  
       END IF                                                                   
c       open file for point spread function                                     
        IF (iPSF.NE.0.AND.psftype.LT.3) THEN                                    
c         wavelength dependent PSF are also printed out                         
          CALL Attach(RootName,length,'.psf',fname)                             
          open(23,file=fname,STATUS='UNKNOWN')                                  
          header1 = '    Offset'                                                
          write(23,'(a10,20f10.2)')header1,(lambdaOut(i),i=1,NlambdaOut)        
        END IF                                                                  
c       Spectra for all models in one file '*.stb' if flag=1                    
        IF(iA.eq.1) THEN                                                        
         CALL Attach(RootName,length,'.stb',fname)                              
         open(15,file=fname,STATUS='UNKNOWN')                                   
        END IF                                                                  
c       All radial profiles in one file '*.rtb' if flag=1                       
        IF(iB.eq.1) THEN                                                        
         CALL Attach(RootName,length,'.rtb',fname)                              
         open(16,file=fname,STATUS='UNKNOWN')                                   
        END IF                                                                  
c       All imaging files in  '*.itb' if flag=1                                 
        IF (abs(iC).eq.1) THEN                                                  
         CALL Attach(RootName,length,'.itb',fname)                              
         open(17,file=fname,STATUS='UNKNOWN')                                   
        END IF                                                                  
        IF(iX.eq.1) THEN                                                        
         CALL Attach(RootName,length,'.mtb',fname)                              
         open(18,file=fname,STATUS='UNKNOWN')                                   
        END IF                                                                  
        IF(iInn.eq.1) THEN                                                      
         CALL Attach(RootName,length,'.err',fname)                              
         open(38,file=fname,STATUS='UNKNOWN')                                   
        END IF                                                                  
c     end if for model=1                                                        
      END IF                                                                    
c -------------------------------------------------------------                 
c     the following files are open for EVERY model                              
                                                                                
c     (the headers for .s## and .r## files are moved to PrOut, MN)              
c     open the spectrum file RootName.s##  (## = model number)                  
      IF(iA.GT.1) THEN                                                          
       write(ch5,'(a2,I3.3)') '.s', model                                       
       CALL Attach(RootName,length,ch5,fname)                                   
       open(15,file=fname,STATUS='UNKNOWN')                                     
       IF(denstyp.EQ.0) THEN                                                    
         IF(iA.eq.3) THEN                                                       
           write(ch5,'(a2,I3.3)') '.z', model                                   
           CALL Attach(RootName,length,ch5,fname)                               
           open(25,file=fname,STATUS='UNKNOWN')                                 
         END IF                                                                 
       END IF                                                                   
      END IF                                                                    
c     open the file RootName.r## (y-dependent quantities)                       
      IF(iB.GT.1) THEN                                                          
        write(ch5,'(a2,I3.3)') '.r', model                                      
        CALL Attach(RootName,length,ch5,fname)                                  
        open(16,file=fname,STATUS='UNKNOWN')                                    
      END IF                                                                    
c     open the file RootName.i## (surface brightness)                           
      IF(abs(iC).GT.1) THEN                                                     
        write(ch5,'(a2,I3.3)') '.i', model                                      
        CALL Attach(RootName,length,ch5,fname)                                  
        open(17,file=fname,STATUS='UNKNOWN')                                    
      END IF                                                                    
c     open the file RootName.c## (convolved images) for flag iC<0               
c      if iC=-3 - in a separate file fname.c##                                  
      IF(iC.EQ.-3.AND.iPSF.GT.0) THEN                                           
        write(ch5,'(a2,I3.3)') '.c', model                                      
        CALL Attach(RootName,length,ch5,fname)                                  
        open(21,file=fname,STATUS='UNKNOWN')                                    
      END IF                                                                    
c     open the file RootName.v## (visibility curves)                            
      IF((abs(iC).EQ.3).AND.(iV.GT.0)) THEN                                     
        write(ch5,'(a2,I3.3)') '.v', model                                      
        CALL Attach(RootName,length,ch5,fname)                                  
        open(22,file=fname,STATUS='UNKNOWN')                                    
      END IF                                                                    
c     open the output file RootName.m##                                         
      IF (iX.GT.1) THEN                                                         
        write(ch5,'(a2,I3.3)') '.m', model                                      
        CALL Attach(RootName,length,ch5,fname)                                  
        open(18,file=fname,STATUS='UNKNOWN')                                    
      END IF                                                                    

      RETURN                                                                    
      END                                                                       
                                                                                

      SUBROUTINE PROUT(model)                                                   
c     =======================================================================       
c     This subroutine prints the results out.        [ZI,Feb'96; MN,Mar'99]         
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
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         
      CHARACTER*72 SC21,SC31,SC43,SB98,SB11,Sbet1,Sbet2,STemp,Serr,             
     &             hdint, hdcon, hdvis, s1, su1, s2, su2, Tstr*10               
      CHARACTER*132 hdsp1,hdsp2,hdrslb1,hdrslb2,hdrsph1,hdrsph2,hdrdyn          
      INTEGER iY, iL, i, model, j, unt                                          
      DOUBLE PRECISION PSFN, PSFfunc(20,1000), Elems(25,200),  ETA,             
     &        faux(npL), fsRbol(npY), tht1, xs, xds, xde, g, res, fnorm         
c ----------------------------------------------------------------------        
c     SmC(1..5, model) are found in ANALYSIS                                    
      SmC(6,model) = Td(1,nY)                                                   
c     SmC(7,model) and SmC(8,model) are Teff(L) and Teff(R), respectively,      
c     the Teff of the slab illumination sources, found from F1 in Analysis      
      tht1 = 412.6/(dsqrt(F1))                                                  
c     colors                                                                    
      CALL getFS(SpecChar(9,model),2,0,SC21)                                    
      CALL getFS(SpecChar(10,model),2,0,SC31)                                   
      CALL getFS(SpecChar(11,model),2,0,SC43)                                   
c     error in %                                                                
      IF (SmC(5,model).LT.0.1) THEN                                             
        CALL getFS(SmC(5,model)*100,0,0,Serr)                                   
      ELSE IF (SmC(5,model).LT.1.0) THEN                                        
        CALL getFS(SmC(5,model)*100,0,1,Serr)                                   
      ELSE                                                                      
        CALL getFS(SmC(5,model)*100,0,2,Serr)                                   
      END IF                                                                    
c     dust temperature at y=Y                                                   
      IF (SmC(6,model).LT.99.5) THEN                                            
        CALL getFS(SmC(6,model),0,0,STemp)                                      
      ELSE                                                                      
        CALL getFS(SmC(6,model),0,1,STemp)                                      
      END IF                                                                    
                                                                                
c --------------  overall parameters to *.OUT file -----------------------      
c     write header to output file *.OUT                                         
      IF (model.EQ.1) THEN                                                      
       write(12,*)'         '                                                   
       write(12,*)' RESULTS:'                                                   
       write(12,*)' --------'                                                   
       IF (denstyp.EQ.0) THEN                                                   
c     ---------- slab output ----------------                                   
        s1=' ###   tau0    Fe1(W/m2)   f1     r1(cm)  Td(K)  Te(L)'             
       su1=' ###     1        2         3       4       5      6  '             
        IF(ksi.GT.0) THEN                                                       
          s2='   Te(R)   err'                                                   
         su2='     7      8 '                                                   
         write(12,'(a54,a14)')s1,s2                                             
         write(12,'(a54,a14)')su1,su2                                           
         su1=' ====================================================='           
         su2='=============='                                                   
         write(12,'(a54,a14)')su1,su2                                           
        ELSE                                                                    
          s2='  err '                                                           
         su2='   7  '                                                           
         write(12,'(a54,a6)')s1,s2                                              
         write(12,'(a54,a6)')su1,su2                                            
         su1=' ====================================================='           
         su2='====='                                                            
         write(12,'(a54,a6)')su1,su2                                            
        END IF                                                                  
       ELSE                                                                     
c     ------------- output for sphere -----------------------------             
        s1=' ###   tau0   F1(W/m2)  r1(cm)    r1/rc   theta1  Td(Y) err'        
       su1=' ###     1        2        3        4        5      6    7 '        
        IF(denstyp.eq.4.OR.RDW) THEN                                            
        s2='   Mdot      Ve       M> '                                          
       su2='     8        9       10 '                                          
         write(12,'(a59,a25)')s1,s2                                             
         write(12,'(a59,a25)')su1,su2                                           
       su1=' =========================================================='        
       su2='============================'                                       
         write(12,'(a59,a28)')su1,su2                                           
        ELSE                                                                    
         write(12,'(a59)')s1                                                    
         write(12,'(a59)')su1                                                   
       su1=' =========================================================='        
         write(12,'(a59)')su1                                                   
        END IF                                                                  
       END IF                                                                   
      END IF                                                                    
c     print output tables                                                       
c     for slab:                                                                 
      IF(denstyp.EQ.0) THEN                                                     
      IF (fmed.LE.accFbol) fmed = 0.                                            
        IF(ksi.NE.0.) THEN                                                      
          CALL Bolom(fsR,fsRbol)                                                
          g = ksi*fsRbol(1) + fmbol(1)                                          
        ELSE                                                                    
          g = fmbol(1)                                                          
        END IF                                                                  
        IF (ksi.GT.0) THEN                                                      
         IF (SmC(6,model).GE.999.5) THEN                                        
          write(12,'(i4,1p,2E9.2,E10.2,E9.2,a5,1p,2E9.2,a4)')                   
     &    model, TAUfid, F1, fmed, Cr1, STemp, SmC(7,model),                    
     &                                         SmC(8,model), Serr               
         ELSE                                                                   
          write(12,'(i4,1p,2E9.2,E10.2,E9.2,a1,a4,1p,2E9.2,a4)')                
     &    model, TAUfid, F1, fmed, Cr1, ' ',STemp, SmC(7,model),                
     &                                             SmC(8,model), Serr           
         END IF                                                                 
        ELSE                                                                    
         IF (SmC(6,model).GE.999.5) THEN                                        
          write(12,'(i4,1p,2E9.2,E10.2,E9.2,a5,1p,E9.2,a3)')                    
     &    model, TAUfid, F1, fmed, Cr1, STemp, SmC(7,model), Serr               
         ELSE                                                                   
          write(12,'(i4,1p,2E9.2,E10.2,E9.2,a1,a4,1p,E9.2,a3)')                 
     &    model, TAUfid, F1, fmed, Cr1,' ', STemp, SmC(7,model), Serr           
         END IF                                                                 
        END IF                                                                  
                                                                                
      ELSE                                                                      
c     for spherical shell                                                       
        IF (SmC(6,model).LT.9.5) THEN                                           
          IF (denstyp.eq.4.OR.RDW) THEN                                         
           write(12,'(i4,1p,5E9.2,a2,a3,a1,a3,1p,3E9.2)')                       
     &     model, TAUfid, F1, Cr1, r1rs, tht1,'  ', STemp,' ',Serr,             
     &                                    CMdot, CVe, CM                        
          ELSE                                                                  
           write(12,'(i4,1p,5E9.2,a2,a4,a3)')                                   
     &     model, TAUfid, F1, Cr1, r1rs, tht1,'  ', STemp, Serr                 
          END IF                                                                
        ELSE                                                                    
         IF (denstyp.eq.4.OR.RDW) THEN                                          
          write(12,'(i4,1p,5E9.2,2(a1,a4),1p,3E9.2)')                           
     &    model, TAUfid, F1, Cr1, r1rs, tht1,' ',STemp,' ', Serr,               
     &                                    CMdot, CVe, CM                        
         ELSE                                                                   
          write(12,'(i4,1p,5E9.2,2(a1,a4))')                                    
     &    model, TAUfid, F1, Cr1, r1rs, tht1,' ',STemp,' ', Serr                
         END IF                                                                 
        END IF                                                                  
                                                                                
       IF (startyp(1).EQ.1.OR.startyp(1).EQ.2) THEN                             
        IF(Tstar.LT.Te_min) THEN                                                
          CALL GetFS(Tstar,0,1,Tstr)                                            
          write(12,'(a50,a5,a5)')                                               
     &' ** WARNING: the input spectrum is a black-body at ',Tstr,' K **'        
          CALL GetFS(Te_min,0,1,Tstr)                                           
          write(12,'(a50,a5,a5)')                                               
     &' *The point-source assumption requires min Teff of ',Tstr,' K **'        
        END IF                                                                  
       END IF                                                                   
c     end if for geometry                                                       
      END IF                                                                    
                                                                                
      CALL getFS(SpecChar(1,model),2,0,SB98)                                    
      CALL getFS(SpecChar(2,model),2,0,SB11)                                    
      CALL getFS(SpecChar(4,model),2,0,Sbet1)                                   
      CALL getFS(SpecChar(5,model),2,0,Sbet2)                                   
c --------------     spectral properties to *.SPP file   -----------------------
      IF (iSPP.NE.0) THEN                                                       
       write(19,'(i3,1p,5E9.2,7a6,E9.2)') model, TAUfid, SmC(1,model),          
     &      (SpecChar(i,model),i=6,8), SC21, SC31, SC43, Sbet1, Sbet2,          
     &       SB98, SB11, SpecChar(3,model)                                      
      END IF                                                                    
                                                                                
c --------------   spectrum to *.s## (old *.Axx) file   ------------------------
      IF (iA.NE.0) THEN                                                         
       IF(denstyp.EQ.0) THEN                                                    
        hdsp1 = '#  lambda     fRight       xAtt       xDs   '                  
       ELSE                                                                     
        hdsp1 = '#  lambda      fTot        xAtt       xDs   '                  
       END IF                                                                   
        hdsp2 = '     xDe       fInp       tauT      albedo  '                  
        unt = 15                                                                
        call line(1,2,unt)                                                      
        IF(denstyp.EQ.0) THEN                                                   
          write(unt,'(a8,i4,a36)')                                              
     &         '#  model',model,'   SPECTRUM from the right slab side '         
        ELSE                                                                    
          write(unt,'(a8,i4,a12)') '#  model',model,'   SPECTRUM '              
        END IF                                                                  
        call line(1,1,unt)                                                      
        DO iL = 1, nL                                                           
c         The right-side spectra are: ftot + ksi*fsR = fsL + fds + fde          
          ftot(iL,nYok) = ftot(iL,nYok) + ksi*fsR(iL,nYok)                      
          faux(iL) = ftot(iL,nYok)/lambda(iL)                                   
        END DO                                                                  
        CALL Simpson(npL,1,nL,lambda,faux,res)                                  
c       Normalization factor for slab s-spectra                                 
        fnorm = res	                                                            
        DO iL = 1, nL                                                           
           IF (ftot(iL,nYok).NE.0.0) THEN                                       
             xs = fsL(iL,nYok)/ftot(iL,nYok)                                    
             xds = fds(iL,nYok)/ftot(iL,nYok)                                   
             xde = fde(iL,nYok)/ftot(iL,nYok)                                   
           ELSE                                                                 
             xs = 0.0                                                           
             xds = 0.0                                                          
             xde = 0.0                                                          
           END IF                                                               
c          change introduced for version 11a - no need to print negligible value
           IF (ftot(iL,nYok).LT.1.0E-20) ftot(iL,nYok) = 0.0                    
           IF (xs.LT.1.0E-20) xs = 0.0                                          
           IF (xds.LT.1.0E-20) xds = 0.0                                        
           IF (xde.LT.1.0E-20) xde = 0.0                                        
           IF (fsL(iL,1).LT.1.0E-20) fsL(iL,1) = 0.0                            
c          for slab rescale fTot with the bolom flux                            
           IF (denstyp.eq.0) ftot(iL,nYok) = ftot(iL,nYok)/fnorm                
           Elems(1,iL) = lambda(iL)                                             
           Elems(2,iL) = fTot(iL,nYok)                                          
           Elems(3,iL) = xs                                                     
           Elems(4,iL) = xds                                                    
           Elems(5,iL) = xde                                                    
           Elems(6,iL) = fsL(iL,1)                                              
           Elems(7,iL) = TAUtot(iL)                                             
           Elems(8,iL) = omega(iL,1)                                            
        END DO                                                                  
c     ------ Tabulate the spectra in the desired form ----------                
        write(unt,'(2(a45))') hdsp1,hdsp2                                       
        CALL MakeTable(Elems,nL,8,unt)                                          
c       spectra from the illuminated slab side (file *.z##)                     
        IF (denstyp.EQ.0) THEN                                                  
           DO iL = 1, nL                                                        
c           The left-side spectra are: |R*fsR + fm|=|fsL-ftot|                  
            ftot(iL,1) = dabs(fsL(iL,1) - ftot(iL,1))                           
c           normalization of z-spectra for slab                                 
            faux(iL) = ftot(iL,1)/lambda(iL)                                    
           END DO                                                               
           CALL Simpson(npL,1,nL,lambda,faux,res)                               
           fnorm = res	                                                         
           DO iL = 1, nL                                                        
             IF (ftot(iL,1).NE.0.0) THEN                                        
                xs =  ksi*fsR(iL,1)/ftot(iL,1)                                  
                xds = dabs(fds(iL,1)/ftot(iL,1))                                
                xde = dabs(fde(iL,1)/ftot(iL,1))                                
             ELSE                                                               
                xs = 0.0                                                        
                xds = 0.0                                                       
                xde = 0.0                                                       
             END IF                                                             
             IF (xs.LT.1.0E-20) xs =0.0                                         
             IF (xds.LT.1.0E-20) xds =0.0                                       
             IF (xde.LT.1.0E-20) xde =0.0                                       
             IF(fsR(iL,nYok).LT.1.0E-20) fsR(iL,nYok) = 0.0                     
             IF(ftot(iL,1).LT.1.0E-20) ftot(iL,1) = 0.0                         
c            rescale fTot with the bolom flux for z-spectra                     
             ftot(iL,1) = ftot(iL,1)/fnorm                                      
             Elems(1,iL) = lambda(iL)                                           
             Elems(2,iL) = fTot(iL,1)                                           
             Elems(3,iL) = xs                                                   
             Elems(4,iL) = xds                                                  
             Elems(5,iL) = xde                                                  
             Elems(6,iL) = ksi*fsR(iL,nYok)                                     
             Elems(7,iL) = TAUtot(iL)                                           
             Elems(8,iL) = omega(iL,1)                                          
           END DO                                                               
           IF (iA.EQ.3) unt=25                                                  
c          append to the .s## file or write in a separate .z## file (if iA=3)   
           call line(1,1,unt)                                                   
           write(unt,'(a8,i4,a36)')                                             
     &         '#  model',model,'   SPECTRUM from the left slab side'           
           call line(1,1,unt)                                                   
           hdsp1 = '#  lambda      fLeft       xAtt       xDs   '               
           write(unt,'(2(a45))') hdsp1,hdsp2                                    
           CALL MakeTable(Elems,nL,8,unt)                                       
        END IF                                                                  
      END IF                                                                    
                                                                                
c     spectral properties for *.ZPP file in slab case                           
      IF (denstyp.eq.0.AND.iSPP.NE.0) THEN                                      
        CALL getFS(SpecChar(12,model),2,0,SB98)                                 
        CALL getFS(SpecChar(13,model),2,0,SB11)                                 
        CALL getFS(SpecChar(15,model),2,0,Sbet1)                                
        CALL getFS(SpecChar(16,model),2,0,Sbet2)                                
        CALL getFS(SpecChar(20,model),2,0,SC21)                                 
        CALL getFS(SpecChar(21,model),2,0,SC31)                                 
        CALL getFS(SpecChar(22,model),2,0,SC43)                                 
c       write spectral properties to *.ZPP file                                 
         IF (iSPP.EQ.3) THEN                                                    
            write(24,'(i3,1p,5E9.2,7a6,E9.2)')                                  
     &      model, TAUfid, SmC(1,model), (SpecChar(i,model),i=17,19),           
     &      SC21,SC31,SC43, Sbet1,Sbet2, SB98,SB11,SpecChar(14,model)           
         ELSE                                                                   
            write(zline(model),'(i3,1p,5E9.2,7a6,E9.2)')                        
     &      model, TAUfid, SmC(1,model), (SpecChar(i,model),i=17,19),           
     &      SC21,SC31,SC43, Sbet1,Sbet2,SB98,SB11,SpecChar(14,model)            
         END IF                                                                 
      END IF                                                                    
                                                                                
c -----------  y-dependent quantities to *.r## (old *.Bxx) file -------------   
      IF (iB.NE.0) THEN                                                         
       hdrslb1= '#     t        tauF      epsilon       Td '                    
       hdrslb2= '      febol      fRbol      fLbol '                            
       hdrsph1= '#     y         eta         t         tauF '                   
       hdrsph2= '     epsilon      Td         rg '                              
        hdrdyn= '         u        drift'                                       
       unt = 16                                                                 
       call line(1,2,unt)                                                       
       IF(denstyp.EQ.0) THEN                                                    
          write(unt,'(a8,i4,a18)') '#  model',model,'  SPATIAL PROFILES'        
       ELSE                                                                     
          write(unt,'(a8,i4,a18)') '#  model',model,'  RADIAL PROFILES '        
       END IF                                                                   
       call line(1,1,unt)                                                       
c    --------- for slab ---------                                               
       IF (denstyp.EQ.0) THEN                                                   
         DO iY = 1, nYok                                                        
           Elems(1,iY) = tr(iY)                                                 
           Elems(2,iY) = tauF(iY)                                               
           Elems(3,iY) = Eps(iY)                                                
           Elems(4,iY) = Td(1,iY)                                               
           Elems(5,iY) = fsbol(iY)                                              
           Elems(6,iY) = fpbol(iY)                                              
           Elems(7,iY) = fmbol(iY)                                              
         END DO                                                                 
c      check for nonsense values:                                               
         DO i = 1, 8                                                            
          DO iY = 1, nYok                                                       
           IF (Elems(i,iY).NE.Elems(i,iY)) THEN                                 
              Elems(i,iY) = 0.0                                                 
           END IF                                                               
          END DO                                                                
         END DO                                                                 
         write(unt,'(a42,a34)') hdrslb1,hdrslb2                                 
         CALL MakeTable(Elems,nYok,7,unt)                                       
       ELSE                                                                     
c     ------  for spherical shell --------                                      
        DO iY = 1, nYok                                                         
          Elems(1,iY) = Yok(iY)                                                 
          Elems(2,iY) = ETA(Y(iY))                                              
          Elems(3,iY) = tr(iY)                                                  
          Elems(4,iY) = tauF(iY)                                                
          Elems(5,iY) = Eps(iY)                                                 
          Elems(6,iY) = Td(1,iY)                                                
          Elems(7,iY) = rg(1,iY)                                                
        END DO                                                                  
c       check for nonsense values:                                              
        DO i = 1, 7                                                             
         DO iY = 1, nYok                                                        
          IF(Elems(i,iY).NE.Elems(i,iY).OR.Elems(i,iY).LT.1.e-20) THEN          
            Elems(i,iY) = 0.0                                                   
          END IF                                                                
         END DO                                                                 
        END DO                                                                  
c       with dynamics                                                           
        IF (RDW) THEN                                  
         DO iY = 1, nYok                                                        
           Elems(8,iY) = ugas(iY)/ugas(nYok)                                    
           Elems(9,iY) = vrat(1,iY)                                             
         END DO                                                                 
c        check values:                                                          
         DO i = 8, 9                                                            
          DO iY = 1, nYok                                                       
           IF(Elems(i,iY).NE.Elems(i,iY).OR.Elems(i,iY).LT.1.e-20) THEN         
             Elems(i,iY) = 0.0                                                  
           END IF                                                               
          END DO                                                                
         END DO                                                                 
         write(unt,'(a42,a32,a23)') hdrsph1,hdrsph2,hdrdyn                      
         CALL MakeTable(Elems,nYok,9,unt)                                       
        ELSE                                                                    
         write(unt,'(a42,a32)') hdrsph1,hdrsph2                                 
         CALL MakeTable(Elems,nYok,7,unt)                                       
        END IF                                                                  
c      end if for geometry                                                      
       END IF                                                                   
c     end if for the iB (radial) flag                                           
      END IF                                                                    
                                                                                
c --------------   intensities to *.inn (old *.Cxx) file  --------------        
      IF (abs(iC).NE.0) THEN                                                    
        hdint = '#     b          t(b)'                                         
        hdcon = '#   Offset '                                                   
        hdvis = '#     q    '                                                   
        unt = 17                                                                
        CALL LINE(1,2,unt)                                                      
        write(unt,'(a8,i4,a14)') '#  model',model,'   RAW IMAGE  '              
        CALL LINE(1,1,unt)                                                      
        DO i = 1, nPok+2                                                        
          Elems(1,i) = bOut(i)                                                  
          Elems(2,i) = tauZout(i)                                               
          DO j = 1, NLambdaOut                                                  
c           check values:                                                       
            IF(IntOut(j,i).NE.IntOut(j,i).OR.IntOut(j,i).LT.1.e-20) THEN        
               IntOut(j,i) = 0.0                                                
            END IF                                                              
            Elems(j+2,i) = IntOut(j,i)                                          
c           we want intensity in Jy/arcsec2                                     
            Elems(j+2,i) = 7.83 * lambdaOut(j) * F1 * Elems(j+2,i)              
          END DO                                                                
        END DO                                                                  
        write(unt,'(a21,20f11.2)')hdint,(lambdaOut(j),j=1,NlambdaOut)           
        CALL MakeTable(Elems,nPok+2,NlambdaOut+2,unt)                           
      END IF                                                                    
      IF (iC.LT.0) THEN                                                         
c ---------  convolved images either add to .i## file or write in *.c## file -- 
        IF(iC.EQ.-3) unt = 21                                                   
        call line(1,2,unt)                                                      
        write(unt,'(a8,i4,a20)') '#  model',model,'   CONVOLVED IMAGE  '        
        call line(1,1,unt)                                                      
        DO i = 1, Nconv                                                         
         Elems(1,i) = Offset(i)                                                 
         DO j = 1, NLambdaOut                                                   
c        check values:                                                          
         IF(ConvInt(j,i).NE.ConvInt(j,i).OR.ConvInt(j,i).LT.1.e-20) THEN        
            ConvInt(j,i) = 0.0                                                  
         END IF                                                                 
         Elems(j+1,i) = ConvInt(j,i)                                            
         END DO                                                                 
        END DO                                                                  
        write(unt,'(a11,20f11.2)')hdcon,(lambdaOut(i),i=1,NlambdaOut)           
        CALL MakeTable(Elems,Nconv,NLambdaOut+1,unt)                            
        IF (psftype.LT.3.AND.model.EQ.1) THEN                                   
c         Wavelength dependent PSFs, print them separately in *.psf             
c         first generate wavelength dependent PSFs                              
          DO j = 1, NlambdaOut                                                  
            iLambda = j                                                         
            DO i = 1, Nconv                                                     
              PSFfunc(j,i) = PSFN(Offset(i))                                    
c             check dynamic range                                               
              CALL CHKRANGE(dynrange,PSFfunc(j,i))                              
            END DO                                                              
          END DO                                                                
c         print them out                                                        
          DO i = 1, Nconv                                                       
           write(23,'(1p,e12.5,20e10.3)')Offset(i),                             
     &                                 (PSFfunc(j,i),j=1,NlambdaOut)            
          END DO                                                                
        END IF                                                                  
      END IF                                                                    
c --------------     visibility curves to *.vnn file    ------------------------
      IF (iV.NE.0) THEN                                                         
        IF(abs(iC).EQ.3) unt = 22                                               
        call line(1,2,unt)                                                      
        write(unt,'(a8,i4,a14)') '#  model',model,'  VISIBILITY  '              
        call line(1,1,unt)                                                      
        DO i = 1, Nvisi                                                         
          Elems(1,i) = qtheta1(i)                                               
          DO j = 1, NLambdaOut                                                  
c          check for nonsense values:                                           
           IF(Visib(j,i).NE.Visib(j,i).OR.Visib(j,i).LT.1.e-20) THEN            
              Visib(j,i) = 0.0                                                  
           END IF                                                               
           Elems(j+1,i) = Visib(j,i)                                            
          END DO                                                                
        END DO                                                                  
        write(unt,'(a11,20f11.2)')hdvis,(lambdaOut(i),i=1,NlambdaOut)           
        CALL MakeTable(Elems,Nvisi,NLambdaOut+1,unt)                            
      END IF                                                                    

      RETURN                                                                    
      END                                                                       

                                                                                
                                                                                

c This function is taken from Moshe Elitzur.           [Z.I., Nov. 1995]        
c =======================================================================       
      DOUBLE PRECISION FUNCTION RDINP(Equal,IUNIT)                              
c =======================================================================       
C     Read lines, up to 232 long, from pre-opened unit IUNIT and extract        
C     all input numbers from them. When EQUAL is set, numeric input data        
C     must be preceded by an equal sign.All non-numeric data and numbers        
C     not preceded by = when EQUAL is on are ignored.RDINP = next number        
C     encountered (after equal sign) and terminated by a nun-numeric            
C     symbol (termination with blank is best). Commas and exponential           
C     notation are allowed.  All text after % is ignored, as in TeX.            
C     Lines with * in the first column are echoed to the output device.         
C     The number is comprised of an actual VALUE, decimal part FRAC and         
C     (integer) exponent PWR.  It has a SIGN, and the exponential power         
C     has sign SIGNEX. Logical flag to the decimal part is set in               
C     DECIMAL. The search is conducted between FIRST, which is                  
C     continuously increased, and LAST.  A new line is read when FIRST          
C     exceeds LAST, and can be forced by calling with -IUNIT.  Actual           
C     extraction of numerical value is done in separate FUNCTION VAL.           
c =======================================================================       
      IMPLICIT None                                                             
      Integer IUnit,Ind,First, Last                                             
      DOUBLE PRECISION Value,Val,Frac,Pwr,Sign,Signex                           
      CHARACTER Card*(232),CR,prev,Term,Next                                    
      Logical Equal,digit,minus,sgn,dot,E,decimal                               
      Save Card,First,Last                                                      
      DATA First/1/, Last/0/                                                    

C FUNCTION STATEMENTS                                                           
      digit(CR) = CR.GE.'0' .AND. CR.LE.'9'                                     
      minus(CR) = CR.EQ.'-'                                                     
      sgn(CR)   = CR.EQ.'+' .OR. CR.EQ.'-'                                      
      dot(CR)   = CR.EQ.'.'                                                     
      E(CR)     = CR.EQ.'E' .OR. CR.EQ.'e'                                      
C                                                                               
      IF (IUnit.lt.0) Then                                                      
         First = Last + 1                                                       
         IUnit = -IUnit                                                         
      END IF                                                                    
c     Start the search for the next number:                                     
  1   RDINP  = 0.                                                               
      VALUE  = 0.                                                               
      FRAC   = 0.                                                               
      PWR    = 0.                                                               
      SIGN   = 1.                                                               
      SIGNEX = 1.                                                               
      Decimal = .False.                                                         
      If (first.gt.last) then                                                   
c     Time to get a new line                                                    
         READ (IUNIT, '(A)' , END = 99) Card                                    
         first = 1                                                              
         last = len(Card)                                                       
c        Find start of trailing junk:                                           
         DO WHILE (Card(last:last).LE.' ')                                      
          last = last - 1                                                       
          if (last.lt.first) goto 1                                             
         END DO	                                                                
         IF (Card(first:first).EQ.'*') WRITE (12,'(A)') Card(1:last)            
         ind = Index(Card,'%')                                                  
         if (ind.gt.0) last = ind - 1                                           
      End If                                                                    
                                                                                
c     Get past the next '=' when the EQUAL flag is set                          
      If (Equal) then                                                           
        DO WHILE (Card(first:first).NE.'=')                                     
          first = first + 1                                                     
          IF (first.GT.last) goto 1                                             
        END DO                                                                  
      End If                                                                    
c     OK, start searching for the next digit:                                   
      Do While (.not.digit(Card(first:first)))                       		         
	    first = first + 1                                                          
		if (first.gt.last) goto 1                                                     
	End Do                                                                         
c     Check if it is a negative or decimal number                               
      If (first.gt.1) then                                                      
         prev = Card(first-1:first-1)                                           
         if (minus(prev)) sign = -1.                                            
         if (dot(prev)) then                                                    
           decimal = .True.                                                     
           if (first.gt.2 .and. minus(Card(first-2:first-2))) sign = -1.        
         end if                                                                 
      End If                                                                    
c     Extract the numerical value                                               
      IF (.not.Decimal) Then                                                    
         Value = VAL(Card,first,last,decimal,Term)                              
c        Check for a possible decimal part.  Termination with '.E' or           
c        '.e' is equivalent to the same ending without '.'                      
         if (first.lt.last.and.dot(Term)) then                                  
            first = first + 1                                                   
            next = Card(first:first)                                            
            if (digit(next)) decimal = .true.                                   
            if (E(next)) Term = 'E'                                             
         end if                                                                 
      END IF                                                                    
c     Extract the decimal fraction, when it exists                              
      IF (Decimal) Frac = VAL(Card,first,last,decimal,Term)                     
c     An exponential may exist if any part terminated with 'E' or 'e'           
      IF (first.lt.last.and.E(term)) then                                       
         first = first + 1                                                      
         next = Card(first:first)                                               
         if (first.lt.last.and.sgn(next))then                                   
            first = first + 1                                                   
            if (minus(next)) Signex = -1.                                       
         end if                                                                 
         decimal = .False.                                                      
         PWR = VAL(Card,first,last,decimal,Term)                                
      END IF                                                                    
c     Finally, put the number together                                          
      RDINP = Sign*(Value + Frac)*10**(Signex*PWR)                              
      Return                                                                    
                                                                                
99    WRITE (12,'(3(1x,a,/))')                                                  
     *' ****TERMINATED. EOF reached by RDINP while looking for input. ',        
     *' *** Last line read:',Card                                               

      RETURN                                                                    
      END                                                                       

                                                                                
                                                                                
                                                                                

c This function is taken from Moshe Elitzur            [Z.I., Nov. 1995]        
c =======================================================================       
      DOUBLE PRECISION FUNCTION VAL(Card,first,last,decimal,Term)               
c     Extract numerical value from CARD, begining at position FIRST up          
c     to the first non-digit encountered.  The terminating character is         
c     returned in TERM, its position in FIRST. Commas are treated as            
c     simple separators.                                                        
c =======================================================================       
      IMPLICIT None                                                             
      Character Card*(*), Term, CH                                              
      Logical Decimal, digit                                                    
      Integer IVAL, first, last, first0                                         
      Real pwr                                                                  

C     FUNCTION STATEMENTS                                                       
      IVAL(CH) = ICHAR(CH) - ICHAR('0')                                         
      digit(CH) = CH.GE.'0' .AND. CH.LE.'9'                                     
c                                                                               
      VAL = 0.                                                                  
      pwr = 1.                                                                  
      first0 = first                                                            
      DO 10 first = first0, last                                                
         Term = Card(first:first)                                               
         if (Term.eq.',') goto 10                                               
         if (.not.digit(Term)) Return                                           
         If (decimal) then                                                      
            pwr = pwr*0.1                                                       
            Val = Val + pwr*IVAL(Term)                                          
         Else                                                                   
            Val = 10.*Val + IVAL(Term)                                          
         End If                                                                 
  10  CONTINUE                                                                  
      Term = ' '                                                                

      RETURN                                                                    
      END                                                                       
