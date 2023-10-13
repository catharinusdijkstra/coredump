                                                                                
c**********************************************************************         
c     This is the block data with optical constants for the supported           
c     grain types. It has to be at the beginning of the program.                
c     The data was compiled from different sources by Z. Ivezic (1996).         
c                                                          [MN, Apr'98]         
c =====================================================================         
      BLOCK DATA                                                                
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
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
      DATA lam_nk/1.00E-02,2.00E-02,3.00E-02,4.00E-02,5.00E-02,                 
     &   6.00E-02,8.00E-02,1.00E-01,1.20E-01,1.50E-01,2.00E-01,2.50E-01,        
     &   3.00E-01,3.60E-01,4.40E-01,5.50E-01,7.00E-01,8.50E-01,1.00E+00,        
     &   1.15E+00,1.30E+00,1.70E+00,2.00E+00,2.20E+00,2.70E+00,3.00E+00,        
     &   3.50E+00,4.00E+00,4.50E+00,5.00E+00,5.50E+00,6.00E+00,6.50E+00,        
     &   7.00E+00,7.50E+00,8.00E+00,8.50E+00,9.00E+00,9.40E+00,9.55E+00,        
     &   9.70E+00,9.85E+00,1.00E+01,1.05E+01,1.10E+01,1.13E+01,1.16E+01,        
     &   1.20E+01,1.25E+01,1.30E+01,1.35E+01,1.40E+01,1.45E+01,1.50E+01,        
     &   1.60E+01,1.70E+01,1.80E+01,1.90E+01,2.00E+01,2.20E+01,2.30E+01,        
     &   2.40E+01,2.50E+01,2.60E+01,2.70E+01,2.80E+01,3.00E+01,3.50E+01,        
     &   4.00E+01,4.50E+01,5.00E+01,5.50E+01,6.00E+01,6.50E+01,7.00E+01,        
     &   7.50E+01,8.00E+01,8.50E+01,9.00E+01,9.50E+01,1.00E+02,1.05E+02,        
     &   1.10E+02,1.20E+02,1.30E+02,1.40E+02,1.50E+02,2.00E+02,2.50E+02,        
     &   3.00E+02,4.00E+02,5.00E+02,7.00E+02,1.30E+03,4.00E+03,1.30E+04,        
     &   2.00E+04,3.60E+04/                                                     
      DATA n_sil_ow/1.81E+00,1.81E+00,1.81E+00,1.81E+00,1.81E+00,               
     &   1.81E+00,1.81E+00,1.81E+00,1.81E+00,1.81E+00,1.81E+00,1.81E+00,        
     &   1.81E+00,1.81E+00,1.82E+00,1.84E+00,1.85E+00,1.85E+00,1.85E+00,        
     &   1.85E+00,1.85E+00,1.87E+00,1.88E+00,1.88E+00,1.89E+00,1.89E+00,        
     &   1.88E+00,1.87E+00,1.87E+00,1.85E+00,1.83E+00,1.81E+00,1.79E+00,        
     &   1.75E+00,1.71E+00,1.64E+00,1.55E+00,1.45E+00,1.40E+00,1.43E+00,        
     &   1.44E+00,1.46E+00,1.51E+00,1.73E+00,1.91E+00,2.00E+00,2.11E+00,        
     &   2.21E+00,2.21E+00,2.16E+00,2.09E+00,2.06E+00,2.03E+00,1.98E+00,        
     &   1.91E+00,1.90E+00,1.96E+00,2.04E+00,2.12E+00,2.27E+00,2.34E+00,        
     &   2.36E+00,2.38E+00,2.40E+00,2.42E+00,2.43E+00,2.45E+00,2.60E+00,        
     &   2.66E+00,2.70E+00,2.74E+00,2.76E+00,2.77E+00,2.79E+00,2.80E+00,        
     &   2.81E+00,2.82E+00,2.83E+00,2.85E+00,2.86E+00,2.87E+00,2.87E+00,        
     &   2.87E+00,2.87E+00,2.88E+00,2.88E+00,2.88E+00,2.89E+00,2.90E+00,        
     &   2.90E+00,2.90E+00,2.90E+00,2.91E+00,2.91E+00,2.91E+00,2.91E+00,        
     &   2.91E+00,2.91E+00/                                                     
      DATA k_sil_ow/9.99E-02,9.99E-02,9.99E-02,9.99E-02,9.99E-02,               
     &   9.99E-02,9.99E-02,9.99E-02,9.99E-02,9.99E-02,9.99E-02,9.99E-02,        
     &   9.99E-02,9.99E-02,9.08E-02,7.17E-02,6.03E-02,5.64E-02,5.46E-02,        
     &   6.12E-02,6.75E-02,7.23E-02,7.33E-02,6.96E-02,6.19E-02,5.88E-02,        
     &   5.48E-02,5.24E-02,5.07E-02,5.14E-02,5.24E-02,5.50E-02,5.79E-02,        
     &   6.36E-02,7.07E-02,9.03E-02,1.31E-01,2.63E-01,4.31E-01,4.96E-01,        
     &   5.60E-01,6.30E-01,7.12E-01,8.10E-01,7.94E-01,7.86E-01,7.68E-01,        
     &   6.13E-01,4.93E-01,3.98E-01,3.74E-01,3.69E-01,3.66E-01,3.72E-01,        
     &   4.39E-01,5.63E-01,6.68E-01,7.23E-01,7.51E-01,7.55E-01,7.34E-01,        
     &   7.10E-01,6.85E-01,6.61E-01,6.37E-01,6.44E-01,6.58E-01,6.28E-01,        
     &   5.78E-01,5.23E-01,4.69E-01,4.34E-01,4.13E-01,3.92E-01,3.70E-01,        
     &   3.49E-01,3.28E-01,3.06E-01,2.85E-01,2.64E-01,2.42E-01,2.37E-01,        
     &   2.32E-01,2.22E-01,2.12E-01,2.02E-01,1.92E-01,1.42E-01,1.00E-01,        
     &   9.09E-02,7.24E-02,5.39E-02,3.79E-02,2.15E-02,7.26E-03,2.39E-03,        
     &   2.39E-03,2.39E-03/                                                     
      DATA n_sil_oc/1.77E+00,1.77E+00,1.77E+00,1.77E+00,1.77E+00,               
     &   1.77E+00,1.77E+00,1.77E+00,1.77E+00,1.77E+00,1.77E+00,1.77E+00,        
     &   1.77E+00,1.77E+00,1.79E+00,1.82E+00,1.83E+00,1.82E+00,1.81E+00,        
     &   1.81E+00,1.81E+00,1.84E+00,1.85E+00,1.86E+00,1.87E+00,1.88E+00,        
     &   1.88E+00,1.87E+00,1.86E+00,1.85E+00,1.83E+00,1.80E+00,1.78E+00,        
     &   1.74E+00,1.69E+00,1.62E+00,1.51E+00,1.39E+00,1.34E+00,1.36E+00,        
     &   1.37E+00,1.39E+00,1.44E+00,1.67E+00,1.86E+00,1.96E+00,2.08E+00,        
     &   2.22E+00,2.23E+00,2.17E+00,2.09E+00,2.04E+00,2.01E+00,1.92E+00,        
     &   1.78E+00,1.73E+00,1.80E+00,1.98E+00,2.17E+00,2.40E+00,2.48E+00,        
     &   2.53E+00,2.58E+00,2.63E+00,2.69E+00,2.70E+00,2.74E+00,2.90E+00,        
     &   2.95E+00,2.98E+00,3.01E+00,3.02E+00,3.03E+00,3.03E+00,3.03E+00,        
     &   3.04E+00,3.04E+00,3.05E+00,3.05E+00,3.05E+00,3.06E+00,3.06E+00,        
     &   3.06E+00,3.06E+00,3.06E+00,3.06E+00,3.06E+00,3.07E+00,3.07E+00,        
     &   3.07E+00,3.07E+00,3.08E+00,3.08E+00,3.08E+00,3.08E+00,3.08E+00,        
     &   3.08E+00,3.08E+00/                                                     
      DATA k_sil_oc/8.95E-02,8.95E-02,8.95E-02,8.95E-02,8.95E-02,               
     &   8.95E-02,8.95E-02,8.95E-02,8.95E-02,8.95E-02,8.95E-02,8.95E-02,        
     &   8.95E-02,8.95E-02,8.23E-02,6.64E-02,5.09E-02,4.46E-02,4.69E-02,        
     &   6.37E-02,7.89E-02,9.62E-02,1.03E-01,1.00E-01,9.27E-02,8.97E-02,        
     &   8.50E-02,8.14E-02,7.85E-02,7.80E-02,7.77E-02,8.06E-02,8.40E-02,        
     &   8.91E-02,9.51E-02,1.10E-01,1.48E-01,2.95E-01,4.78E-01,5.49E-01,        
     &   6.19E-01,6.96E-01,7.87E-01,9.11E-01,9.10E-01,9.14E-01,9.07E-01,        
     &   7.40E-01,5.94E-01,4.74E-01,4.39E-01,4.28E-01,4.20E-01,4.35E-01,        
     &   5.35E-01,7.36E-01,9.52E-01,1.10E+00,1.13E+00,1.07E+00,1.04E+00,        
     &   1.01E+00,0.97E+00,0.94E+00,8.98E-01,8.77E-01,8.42E-01,7.52E-01,        
     &   6.60E-01,5.68E-01,4.75E-01,4.25E-01,4.04E-01,3.82E-01,3.60E-01,        
     &   3.39E-01,3.17E-01,2.96E-01,2.74E-01,2.52E-01,2.31E-01,2.26E-01,        
     &   2.21E-01,2.12E-01,2.02E-01,1.92E-01,1.83E-01,1.35E-01,9.44E-02,        
     &   8.56E-02,6.82E-02,5.07E-02,3.57E-02,2.02E-02,6.82E-03,2.29E-03,        
     &   2.29E-03,2.29E-03/                                                     
      DATA n_sil_dl/8.66E-01,8.66E-01,8.66E-01,8.25E-01,7.81E-01,               
     &   8.90E-01,1.29E+00,1.60E+00,1.87E+00,2.26E+00,1.93E+00,1.80E+00,        
     &   1.76E+00,1.74E+00,1.73E+00,1.72E+00,1.72E+00,1.71E+00,1.71E+00,        
     &   1.71E+00,1.71E+00,1.71E+00,1.71E+00,1.71E+00,1.70E+00,1.70E+00,        
     &   1.69E+00,1.68E+00,1.66E+00,1.64E+00,1.60E+00,1.57E+00,1.52E+00,        
     &   1.47E+00,1.34E+00,1.21E+00,1.16E+00,1.11E+00,1.20E+00,1.24E+00,        
     &   1.30E+00,1.35E+00,1.39E+00,1.57E+00,1.75E+00,1.82E+00,1.90E+00,        
     &   2.00E+00,2.04E+00,2.09E+00,2.04E+00,2.00E+00,1.91E+00,1.82E+00,        
     &   1.73E+00,1.69E+00,1.71E+00,1.78E+00,1.89E+00,2.04E+00,2.11E+00,        
     &   2.18E+00,2.26E+00,2.30E+00,2.35E+00,2.40E+00,2.50E+00,2.63E+00,        
     &   2.76E+00,2.89E+00,3.02E+00,3.07E+00,3.13E+00,3.18E+00,3.23E+00,        
     &   3.25E+00,3.27E+00,3.29E+00,3.31E+00,3.33E+00,3.35E+00,3.35E+00,        
     &   3.36E+00,3.37E+00,3.38E+00,3.39E+00,3.40E+00,3.41E+00,3.42E+00,        
     &   3.42E+00,3.43E+00,3.43E+00,3.43E+00,3.43E+00,3.43E+00,3.43E+00,        
     &   3.43E+00,3.43E+00/                                                     
      DATA k_sil_dl/1.39E-01,1.39E-01,1.39E-01,2.55E-01,4.29E-01,               
     &   6.73E-01,8.78E-01,9.26E-01,7.22E-01,5.30E-01,5.32E-02,2.77E-02,        
     &   2.84E-02,2.88E-02,2.91E-02,2.94E-02,2.97E-02,3.00E-02,3.03E-02,        
     &   3.06E-02,3.09E-02,3.21E-02,3.31E-02,3.39E-02,3.60E-02,3.72E-02,        
     &   3.94E-02,4.11E-02,4.25E-02,4.40E-02,4.72E-02,5.05E-02,5.36E-02,        
     &   5.66E-02,1.14E-01,1.71E-01,3.68E-01,5.66E-01,7.65E-01,8.29E-01,        
     &   8.71E-01,8.97E-01,9.24E-01,9.63E-01,1.00E+00,9.65E-01,9.28E-01,        
     &   8.78E-01,7.71E-01,6.63E-01,5.70E-01,4.77E-01,4.83E-01,4.88E-01,        
     &   5.83E-01,7.14E-01,8.55E-01,9.76E-01,1.05E+00,1.08E+00,1.10E+00,        
     &   1.11E+00,1.13E+00,1.12E+00,1.12E+00,1.12E+00,1.11E+00,1.06E+00,        
     &   1.01E+00,9.66E-01,9.19E-01,8.66E-01,8.13E-01,7.59E-01,7.06E-01,        
     &   6.72E-01,6.38E-01,6.03E-01,5.69E-01,5.35E-01,5.00E-01,4.83E-01,        
     &   4.67E-01,4.33E-01,3.99E-01,3.65E-01,3.32E-01,2.48E-01,2.06E-01,        
     &   1.65E-01,1.32E-01,9.87E-02,7.04E-02,3.83E-02,2.46E-02,2.46E-02,        
     &   2.46E-02,2.46E-02/                                                     
      DATA n_amc_hn/8.40E-01,8.40E-01,8.40E-01,8.40E-01,7.40E-01,               
     &   6.90E-01,9.30E-01,1.53E+00,1.74E+00,1.55E+00,1.22E+00,1.40E+00,        
     &   1.60E+00,1.71E+00,1.78E+00,1.85E+00,1.94E+00,2.01E+00,2.11E+00,        
     &   2.21E+00,2.31E+00,2.50E+00,2.63E+00,2.70E+00,2.82E+00,2.86E+00,        
     &   2.95E+00,3.03E+00,3.04E+00,3.04E+00,3.09E+00,3.15E+00,3.16E+00,        
     &   3.16E+00,3.25E+00,3.35E+00,3.38E+00,3.42E+00,3.47E+00,3.49E+00,        
     &   3.51E+00,3.53E+00,3.55E+00,3.59E+00,3.63E+00,3.65E+00,3.67E+00,        
     &   3.70E+00,3.73E+00,3.75E+00,3.79E+00,3.84E+00,3.86E+00,3.88E+00,        
     &   3.96E+00,3.99E+00,4.10E+00,4.14E+00,4.18E+00,4.26E+00,4.31E+00,        
     &   4.34E+00,4.36E+00,4.42E+00,4.47E+00,4.50E+00,4.57E+00,4.77E+00,        
     &   4.94E+00,5.11E+00,5.25E+00,5.32E+00,5.44E+00,5.49E+00,5.67E+00,        
     &   5.71E+00,5.85E+00,5.90E+00,5.99E+00,5.94E+00,6.32E+00,6.70E+00,        
     &   6.79E+00,6.99E+00,7.17E+00,7.37E+00,7.55E+00,8.50E+00,9.32E+00,        
     &   1.01E+01,1.15E+01,1.25E+01,1.41E+01,1.65E+01,1.65E+01,1.65E+01,        
     &   1.65E+01,1.65E+01/                                                     
      DATA k_amc_hn/1.08E-01,1.08E-01,1.08E-01,1.08E-01,1.77E-01,               
     &   3.80E-01,9.00E-01,8.40E-01,5.60E-01,1.77E-01,3.21E-01,7.40E-01,        
     &   7.20E-01,6.86E-01,6.70E-01,6.95E-01,7.70E-01,8.25E-01,9.00E-01,        
     &   9.38E-01,9.60E-01,9.95E-01,1.02E+00,1.01E+00,9.96E-01,9.90E-01,        
     &   1.01E+00,1.04E+00,1.04E+00,1.03E+00,1.09E+00,1.15E+00,1.20E+00,        
     &   1.25E+00,1.34E+00,1.42E+00,1.44E+00,1.47E+00,1.50E+00,1.51E+00,        
     &   1.52E+00,1.53E+00,1.54E+00,1.57E+00,1.60E+00,1.61E+00,1.62E+00,        
     &   1.63E+00,1.66E+00,1.69E+00,1.71E+00,1.74E+00,1.76E+00,1.78E+00,        
     &   1.83E+00,1.89E+00,1.95E+00,1.95E+00,1.98E+00,2.05E+00,2.08E+00,        
     &   2.13E+00,2.19E+00,2.23E+00,2.27E+00,2.32E+00,2.40E+00,2.60E+00,        
     &   2.77E+00,2.92E+00,3.00E+00,3.18E+00,3.31E+00,3.50E+00,3.70E+00,        
     &   3.80E+00,4.00E+00,4.15E+00,4.30E+00,4.48E+00,4.59E+00,4.70E+00,        
     &   4.79E+00,4.97E+00,5.15E+00,5.33E+00,5.51E+00,6.41E+00,7.17E+00,        
     &   7.92E+00,9.43E+00,1.04E+01,1.25E+01,1.41E+01,1.41E+01,1.41E+01,        
     &   1.41E+01,1.41E+01/                                                     
      DATA n_sic_pg/7.06E-01,7.06E-01,7.06E-01,7.06E-01,7.06E-01,               
     &   7.06E-01,7.06E-01,7.06E-01,9.69E-01,2.07E+00,4.29E+00,5.16E+00,        
     &   4.75E+00,3.42E+00,2.59E+00,2.55E+00,2.51E+00,2.50E+00,2.47E+00,        
     &   2.49E+00,2.50E+00,2.51E+00,2.55E+00,2.58E+00,2.65E+00,2.69E+00,        
     &   2.73E+00,2.76E+00,2.76E+00,2.77E+00,2.74E+00,2.71E+00,2.67E+00,        
     &   2.63E+00,2.57E+00,2.51E+00,2.42E+00,2.33E+00,2.17E+00,2.11E+00,        
     &   2.04E+00,1.96E+00,1.88E+00,1.56E+00,1.25E+00,1.47E+00,1.69E+00,        
     &   1.99E+00,3.07E+00,4.14E+00,4.21E+00,4.28E+00,4.11E+00,3.93E+00,        
     &   3.71E+00,3.59E+00,3.52E+00,3.49E+00,3.46E+00,3.43E+00,3.42E+00,        
     &   3.40E+00,3.39E+00,3.39E+00,3.39E+00,3.38E+00,3.38E+00,3.38E+00,        
     &   3.39E+00,3.39E+00,3.39E+00,3.40E+00,3.41E+00,3.41E+00,3.42E+00,        
     &   3.42E+00,3.43E+00,3.43E+00,3.44E+00,3.44E+00,3.45E+00,3.45E+00,        
     &   3.46E+00,3.46E+00,3.47E+00,3.48E+00,3.49E+00,3.51E+00,3.52E+00,        
     &   3.52E+00,3.52E+00,3.52E+00,3.52E+00,3.52E+00,3.52E+00,3.52E+00,        
     &   3.52E+00,3.52E+00/                                                     
      DATA k_sic_pg/1.53E+00,1.53E+00,1.53E+00,1.53E+00,1.53E+00,               
     &   1.53E+00,1.53E+00,1.53E+00,1.43E+00,1.28E+00,1.06E+00,3.50E-01,        
     &   2.50E-01,1.62E-01,1.04E-01,1.13E-01,1.21E-01,1.30E-01,1.44E-01,        
     &   1.62E-01,1.80E-01,2.17E-01,2.63E-01,2.76E-01,2.77E-01,2.78E-01,        
     &   2.51E-01,2.23E-01,1.88E-01,1.53E-01,1.22E-01,9.08E-02,8.07E-02,        
     &   7.06E-02,6.69E-02,6.31E-02,6.59E-02,6.87E-02,9.61E-02,1.05E-01,        
     &   1.09E-01,1.24E-01,1.38E-01,5.43E-01,9.48E-01,1.40E+00,1.85E+00,        
     &   2.45E+00,2.45E+00,2.45E+00,1.69E+00,9.19E-01,7.00E-01,4.80E-01,        
     &   3.62E-01,2.90E-01,2.80E-01,2.82E-01,2.70E-01,2.57E-01,2.51E-01,        
     &   2.45E-01,2.38E-01,2.35E-01,2.32E-01,2.29E-01,2.23E-01,2.13E-01,        
     &   2.04E-01,1.94E-01,1.84E-01,1.80E-01,1.75E-01,1.71E-01,1.66E-01,        
     &   1.63E-01,1.60E-01,1.57E-01,1.54E-01,1.51E-01,1.48E-01,1.46E-01,        
     &   1.43E-01,1.39E-01,1.34E-01,1.29E-01,1.25E-01,1.01E-01,9.22E-02,        
     &   8.30E-02,8.30E-02,8.30E-02,8.30E-02,8.30E-02,8.30E-02,8.30E-02,        
     &   8.30E-02,8.30E-02/                                                     
      DATA n_gr1_dl/9.81E-01,9.25E-01,8.48E-01,7.72E-01,7.23E-01,               
     &   7.65E-01,9.91E-01,1.39E+00,2.57E+00,1.94E+00,1.61E+00,1.56E+00,        
     &   1.93E+00,2.17E+00,2.35E+00,2.34E+00,2.29E+00,2.25E+00,2.23E+00,        
     &   2.21E+00,2.20E+00,2.19E+00,2.18E+00,2.17E+00,2.17E+00,2.16E+00,        
     &   2.16E+00,2.15E+00,2.15E+00,2.14E+00,2.13E+00,2.13E+00,2.12E+00,        
     &   2.11E+00,2.10E+00,2.09E+00,2.08E+00,2.07E+00,2.06E+00,2.06E+00,        
     &   2.05E+00,2.05E+00,2.04E+00,2.03E+00,2.01E+00,1.98E+00,2.05E+00,        
     &   2.01E+00,1.99E+00,1.97E+00,1.96E+00,1.95E+00,1.93E+00,1.92E+00,        
     &   1.89E+00,1.87E+00,1.84E+00,1.82E+00,1.80E+00,1.76E+00,1.74E+00,        
     &   1.73E+00,1.72E+00,1.72E+00,1.71E+00,1.71E+00,1.71E+00,1.76E+00,        
     &   1.83E+00,1.92E+00,2.02E+00,2.11E+00,2.21E+00,2.30E+00,2.40E+00,        
     &   2.49E+00,2.57E+00,2.66E+00,2.74E+00,2.82E+00,2.90E+00,2.97E+00,        
     &   3.05E+00,3.19E+00,3.33E+00,3.46E+00,3.58E+00,4.15E+00,4.65E+00,        
     &   5.10E+00,5.90E+00,6.59E+00,7.80E+00,9.33E+00,9.33E+00,9.33E+00,        
     &   9.33E+00,9.33E+00/                                                     
      DATA k_gr1_dl/3.08E-03,2.78E-02,9.28E-02,2.09E-01,3.70E-01,               
     &   5.39E-01,9.54E-01,1.65E+00,6.55E-01,1.94E-01,2.75E-01,6.38E-01,        
     &   8.16E-01,6.57E-01,4.79E-01,2.55E-01,1.41E-01,9.41E-02,6.96E-02,        
     &   5.51E-02,4.60E-02,3.33E-02,2.87E-02,2.68E-02,2.42E-02,2.35E-02,        
     &   2.32E-02,2.36E-02,2.44E-02,2.57E-02,2.72E-02,2.91E-02,3.13E-02,        
     &   3.41E-02,3.74E-02,4.14E-02,4.61E-02,5.16E-02,5.65E-02,5.85E-02,        
     &   6.05E-02,6.26E-02,6.47E-02,7.23E-02,8.11E-02,8.83E-02,9.49E-02,        
     &   9.92E-02,1.10E-01,1.22E-01,1.34E-01,1.47E-01,1.61E-01,1.75E-01,        
     &   2.07E-01,2.42E-01,2.79E-01,3.20E-01,3.63E-01,4.55E-01,5.04E-01,        
     &   5.55E-01,6.06E-01,6.59E-01,7.11E-01,7.64E-01,8.68E-01,1.11E+00,        
     &   1.33E+00,1.52E+00,1.69E+00,1.85E+00,1.99E+00,2.11E+00,2.23E+00,        
     &   2.34E+00,2.44E+00,2.54E+00,2.63E+00,2.72E+00,2.81E+00,2.89E+00,        
     &   2.97E+00,3.12E+00,3.27E+00,3.40E+00,3.53E+00,4.12E+00,4.62E+00,        
     &   5.08E+00,5.88E+00,6.58E+00,7.79E+00,9.32E+00,9.32E+00,9.32E+00,        
     &   9.32E+00,9.32E+00/                                                     
      DATA n_gr2_dl/9.80E-01,9.20E-01,8.39E-01,7.64E-01,5.55E-01,               
     &   5.43E-01,1.01E+00,2.08E+00,1.87E+00,1.27E+00,6.88E-01,1.25E+00,        
     &   2.39E+00,2.58E+00,2.66E+00,2.74E+00,2.87E+00,3.03E+00,3.19E+00,        
     &   3.34E+00,3.47E+00,3.76E+00,3.94E+00,4.05E+00,4.29E+00,4.43E+00,        
     &   4.64E+00,4.82E+00,4.99E+00,5.14E+00,5.29E+00,5.43E+00,5.57E+00,        
     &   5.69E+00,5.80E+00,5.90E+00,6.01E+00,6.11E+00,6.20E+00,6.22E+00,        
     &   6.25E+00,6.28E+00,6.31E+00,6.39E+00,6.47E+00,6.52E+00,6.56E+00,        
     &   6.62E+00,6.67E+00,6.74E+00,6.81E+00,6.89E+00,6.96E+00,7.03E+00,        
     &   7.15E+00,7.26E+00,7.36E+00,7.46E+00,7.53E+00,7.83E+00,7.99E+00,        
     &   8.13E+00,8.25E+00,8.47E+00,8.73E+00,9.01E+00,9.60E+00,1.14E+01,        
     &   1.35E+01,1.59E+01,1.88E+01,2.13E+01,2.33E+01,2.48E+01,2.59E+01,        
     &   2.66E+01,2.70E+01,2.70E+01,2.66E+01,2.62E+01,2.56E+01,2.46E+01,        
     &   2.34E+01,2.05E+01,1.75E+01,1.49E+01,1.33E+01,1.33E+01,1.68E+01,        
     &   2.11E+01,3.02E+01,3.89E+01,5.47E+01,7.40E+01,7.40E+01,7.40E+01,        
     &   7.40E+01,7.40E+01/                                                     
      DATA k_gr2_dl/3.09E-03,2.79E-02,9.38E-02,1.36E-01,3.62E-01,               
     &   7.20E-01,1.48E+00,1.24E+00,3.38E-01,2.77E-01,1.09E+00,2.27E+00,        
     &   2.08E+00,1.67E+00,1.56E+00,1.56E+00,1.68E+00,1.82E+00,1.94E+00,        
     &   2.03E+00,2.12E+00,2.33E+00,2.48E+00,2.58E+00,2.82E+00,2.96E+00,        
     &   3.18E+00,3.39E+00,3.60E+00,3.79E+00,4.00E+00,4.19E+00,4.37E+00,        
     &   4.55E+00,4.72E+00,4.91E+00,5.09E+00,5.27E+00,5.41E+00,5.46E+00,        
     &   5.50E+00,5.55E+00,5.60E+00,5.77E+00,5.94E+00,6.04E+00,6.14E+00,        
     &   6.27E+00,6.44E+00,6.63E+00,6.82E+00,7.00E+00,7.17E+00,7.35E+00,        
     &   7.69E+00,8.05E+00,8.42E+00,8.80E+00,9.24E+00,1.01E+01,1.05E+01,        
     &   1.10E+01,1.15E+01,1.20E+01,1.24E+01,1.29E+01,1.38E+01,1.59E+01,        
     &   1.76E+01,1.90E+01,1.95E+01,1.92E+01,1.83E+01,1.73E+01,1.62E+01,        
     &   1.50E+01,1.38E+01,1.26E+01,1.16E+01,1.08E+01,9.98E+00,9.24E+00,        
     &   8.73E+00,8.70E+00,9.94E+00,1.24E+01,1.54E+01,2.94E+01,3.95E+01,        
     &   4.78E+01,6.08E+01,7.09E+01,8.61E+01,1.03E+02,1.03E+02,1.03E+02,        
     &   1.03E+02,1.03E+02/                                                     
                                                                                
      END                                                                       
c =====================================================================         
                                                                                
c **********************************************************************        
      PROGRAM Dusty                                                             
C     ========================================
C     DUSTY Radiative Transfer Code
C     ========================================
C     
C     This program solves the continuum radiative transfer problem for a
C     spherically symmetric envelope or for a plane-parallel slab. All input
C     data are read from files named *.inp, listed in the master input file dusty.inp.
C     For detailed information, refer to the Manual.
C     
C     Authors: [Z.I. and M.N.]
C
C ========================================
      IMPLICIT NONE
      CHARACTER*3 version                                                                                                                          
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
      CHARACTER*230 path, apath                                                 
      CHARACTER*235 nameIn, nameOut, nameQ(npG), nameNK(10)                     
      INTEGER error, nG, model, Nmodel, GridType, io1, Empty, lpath,            
     &        iL, Nrec, Nlam                                                    
C     Nrec is the max number of records for TAUgrid in a file                   
      PARAMETER (Nrec = 1000)                                                   
      DOUBLE PRECISION ETAzp(npP,npY), TAUin(Nrec), TAU1, TAU2, RDINP           
      LOGICAL Equal                                                             
      Equal = .TRUE.                                                            
C ----------------------------------------------------------------------        
C     **************************                                                
C     *** ABOUT THIS VERSION ***                                                
C     **************************                                                
      version= '2.01'                                                            
C     Updated versions of Dusty(2.0) with minor changes and bug fixes
C     start from 2.01.

C     Version (2.0) is the first public release. The code has been              
C     significantly improved in terms of speed and I/O options. All             
C     suggestions of the users of version(1.0) had been taken into              
C     consideration. Finished Oct.'99                                           
                                                                                
C     Version 1.0 is a beta version sent to a few people before the             
C     first public release.  Finished Nov,'96.                                  
C **********************************************************************        
C     *** MAIN ***                                                              
C     first open the file with lambda grid and check that the grid satisfies    
C     certain conditions (the wavelengths are in microns) :                     
      OPEN(4, FILE='lambda_grid.dat', STATUS = 'OLD')                           
      Nlam = RDINP(Equal,4)                                                     
      IF (Nlam.NE.npL) THEN                                                     
        WRITE(*,*)' *************** A BIG ERROR !!! ***************** '          
        WRITE(*,*)'  The number of wavelengths in lambda_grid.dat is  '          
        WRITE(*,*)'  not equal to the specified npL in userpar.inc    '          
        WRITE(*,*)'  Make sure the numbers are the same, recompile    '          
        WRITE(*,*)'  and try again.                                   '          
        WRITE(*,*)' ************************************************* '          
        GOTO 999                                                                 
      END IF                                                                    
C     Initialize lambda array                                                   
      READ(4,*,END=99) (lambda(iL), iL = 1, npL)                                
 99   CLOSE(4)                                                                 
      CALL SORT(lambda,npL)                                                     
C     Check the ends of the lambda grid :                                       
      IF(lambda(1).GT.0.01) THEN                                                
        WRITE(*,*)' *************** WARNING! ********************** '            
        WRITE(*,*)'  The shortest wavelength in lambda_grid.dat has '            
        WRITE(*,*)'  to be 0.01 microns. Correct this and try again!'            
        WRITE(*,*)' *********************************************** '            
        GOTO 999                                                                 
      END IF                                                                    
      IF(lambda(npL).LT.36000.) THEN                                            
        WRITE(*,*)' *************** WARNING! ******************* '               
        WRITE(*,*)'  The longest wavelength in lambda_grid.dat   '               
        WRITE(*,*)'  has to be 36 mm. Correct this and try again!'               
        WRITE(*,*)' ******************************************** '               
        GOTO 999                                                                 
      END IF                                                                    
C     Check the resolution:                                                     
      DO iL = 2, npL                                                            
        IF (lambda(iL)/lambda(iL-1).GT.1.51) THEN                               
          WRITE(*,*)' ***************** WARNING!  *******************'          
          WRITE(*,*)' The ratio of two consecutive wavelengths in the'          
          WRITE(*,*)' grid has to be no bigger than 1.5. You have    '          
          WRITE(*,'(2(a4,1p,e8.2))') '    ',lambda(iL)/lambda(iL-1),            
     &                                ' at ', lambda(iL)                        
          WRITE(*,*)' Correct this and try again!                    '          
          WRITE(*,*)' ***********************************************'          
          GOTO 999                                                              
        END IF                                                                  
      END DO                                                                    
C     open master input file dusty.inp                                          
      OPEN(13,ERR=998,FILE='dusty.inp',STATUS='OLD')                            
      io1 = 0                                                                   
C     loop over input files                                                     
      DO WHILE (io1.GE.0)                                                       
C       read a line from master input file using                                
100     READ(13,'(a)',iostat=io1) apath                                         
        IF(io1.LT.0) THEN                                                       
          STOP                                                                   
        END IF                                                                  
        CALL CLEAN(apath, path, lpath)                                          
C       if not EOF and if line is not empty, or commented, proceed              
        IF (EMPTY(path).NE.1) THEN                                              
C         get input/output file names                                           
          CALL ATTACH(path,lpath,'.inp',nameIn)                                 
          CALL ATTACH(path,lpath,'.out',nameOut)                                
C         read input data                                                       
          CALL INPUT(nameIn,nG,nameOut,nameQ,nameNK,                            
     &               TAU1,TAU2,TAUin,Nrec,GridType,Nmodel,error,version)        
           IF (iVerb.GT.0)                                                      
     &          WRITE(*,'(a24,a80)') ' Working on Input File: ',nameIn          
           IF (iVerb.EQ.2) WRITE(*,*) 'Done with Reading Input'                 
C         if an error reading files go to the next input file                   
C         error=3 means some files are missing                                  
          IF (error.EQ.3) GOTO 100                                              
C         get optical properties                                                
          CALL GETOPTPR(nG,nameQ,nameNK,error)                                  
C         if an error reading files go to the next input file                   
          IF (error.EQ.3) GOTO 100                                              
          IF (iVerb.EQ.2) WRITE(*,*) 'Done with GETOPTPR'                       
C         solve for every model                                                 
          DO model = 1, Nmodel                                                  
           IF (iVerb.GT.0)  WRITE(*,'(a9,i4)') ' model = ',model                
           IF (error.EQ.0) THEN                                                 
C            open output files                                                  
             CALL OPPEN(model,path,lpath)                                       
C            calculate optical depth for current model                          
             CALL GETTAU(model,nG,TAU1,TAU2,TAUin,Nrec,GridType,Nmodel)         
             IF (iVerb.EQ.2)                                                    
     &          WRITE(*,*) 'Done with GETTAU. Going to Solve'                   
C            solve                                                              
             CALL SOLVE(model,nG,error,ETAzp)                                   
             IF (error.EQ.0) THEN                                               
C              calculate spectral characteristics                               
               CALL SPECTRAL(model,denstyp,nL,Lambda)                           
C              write results out                                                
               CALL PROUT(model)                                                
             ELSE                                                               
              GOTO 100                                                          
             END IF                                                             
C            close output files                                                 
             CALL CLLOSE(error,model,Nmodel)                                    
           END IF                                                               
           IF (iVerb.EQ.2) WRITE(*,*) ' ----------- '                           
           IF (iVerb.EQ.2) WRITE(*,*) '    '                                    
C         end of the loop over models                                           
          END DO                                                                
        END IF                                                                  
C     end of the loop over input files                                          
      END DO                                                                    
      IF (iVerb.GT.0) WRITE(*,*) ' End of Input Files '                         
      CLOSE(13)                                                                 
C     end this run                                                              
      GOTO 999                                                                  
C     to execute if the master input file is missing                            
998   WRITE(*,*)' *********** FATAL ERROR IN DUSTY ***********'                 
      WRITE(*,*)' * Master input file dusty.inp is missing!? *'                 
      WRITE(*,*)' ********************************************'                 
C ----------------------------------------------------------------------        
999   STOP                                                                      
      END                                                                       
