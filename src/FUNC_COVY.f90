!file FUNC_COVY.f90
!
!-------------------------------------------------------------------
!Software developed by Alberto Bellin  
!Universita' di Trento,
!Dipartimento di Ingegneria Civile ed Ambientale,
!38123-I Mesiano di Povo, TRENTO
!
!e-mail:   Alberto.Bellin@unitn.it 
!
!Copyright: Alberto Bellin and Yoram Rubin.
!
!
!
!Summary:  This function provides different types of covariance
!          functions, namely: the exponential, the Gaussian (bell shape)
!          the Whittle (see Mizell et al. WRR 18(4), 1053-1067, 1982),
!          the Mizell B (see Mizell et al. WRR 18(4), 1053-1067, 1982). 
!
!
!Package Version: 2.0 April 1997
!
!
!please cite the following paper in papers or reports that use
!the present software:
!
!Bellin A., Y. Rubin, Hydro_gen: A new random field generator
!for correlated properties, Stochastic Hydrology and Hydraulics,
!10(4), 1996.
!
!
!Permission is  granted to anyone to use and modify this packages provided
!that:
!i) the authors are acknowledged by citing the abofe referenced paper;
!ii) the use in any kind of research or job will be cited in the relative
!papers or reports;
!iii) the use of the package is under the user responsability
!NO WARRANTY is given concerning bugs and errors.
!iv) The use or distribution must be free of charge.                 
!v) The package uses the following libraries:
!       a) LINPACK by J. J. Dongarra, J. R. Bunch, C. B. Moler 
!          e G.W. Stewart, for the linear system solution
!       b) BLAS, for linear algebra
!       d) RANLIB by Barry W. Brown and James Lovato,
!          Department of Biomathematics, Box 237
!          the University of Texas, M.D. Anderson Cancer Center
!          1515 Holcombe Boulevard, Huston, TX 77030, for the generation
!          of independent normally distributed random numbers.
!       e) Numerical Recipes by W. H. Press, B. P. Flannery, S. A.
!          Teukolsky, W. T. Vetterling, for the function computing
!          the Bessel Function
!   
!   Copyright conditions of the above referenced libraries are 
!   extended to hydro_gen  
!  
!   Bug reports and hints are welcomed to the following e-mail address:
!   Alberto.Bellin@unitn.it
!------------------------------------------------------------------------------
! 
MODULE MOD_FUNC_COVY
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
INTERFACE
      FUNCTION FUNC_COVY(x,y,c0,omega2,scl,indec)
         IMPLICIT NONE
         REAL*8         ,INTENT(IN) ::x
         REAL*8         ,INTENT(IN) ::y
         REAL*8         ,INTENT(IN) ::c0
         REAL*8         ,INTENT(IN) ::omega2
         REAL*8         ,INTENT(IN) ::scl   
         INTEGER        ,INTENT(IN) ::indec
         !
         REAL*8                     ::func_covy
      ENDFUNCTION
ENDINTERFACE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
ENDMODULE MOD_FUNC_COVY
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
FUNCTION FUNC_COVY(x,y,c0,omega2,scl,indec)
! USE MOD_FUNC_BESSI0
! USE MOD_FUNC_BESSI1
! USE MOD_FUNC_BESSK0
! USE MOD_FUNC_BESSK1
 
 IMPLICIT NONE
 REAL*8         ,INTENT(IN) ::x
 REAL*8         ,INTENT(IN) ::y
 REAL*8         ,INTENT(IN) ::c0
 REAL*8         ,INTENT(IN) ::omega2
 REAL*8         ,INTENT(IN) ::scl         
 INTEGER        ,INTENT(IN) ::indec
 !
 REAL*8                     ::func_covy
 !-----------------------------------------
 REAL*8,parameter  :: pi=datan(1.d0)
 REAL*8            :: r2,r
 
 real*8            :: bessk0,bessk1
 
 !
 !
 !
 
 r=0.0d0
 
 
 r2=x*x+y*y
 !
 IF(r2.gt.0.0)THEN
    r=dsqrt(r2)

 !
 !-----  exponential covariance function -----
 ! 
 
 IF(indec.eq.1) func_covy=dexp(-r) 
 
 
!  !
!  !-----  Gaussian covariance function -----
!  !
  IF(indec.eq.2) func_covy=dexp(-r2)
!  !
  IF(indec.eq.3)THEN
!     !
!     !-----  Whittle covariance function 
!     !       (Mizell et al. WRR 18(4), 1053-1067, 1982) -----
     r=r*pi/(dble(2.0)*scl)
     func_covy=r*bessk1(r)              
  ENDIF
!  !
  IF(indec.eq.4)THEN
!     !
!     !----- Mizell covariance function (Type B) -----
!     !    (Mizell et al. WRR 18(4), 1053-1067, 1982) -----
!     !
     r=dble(3.)*pi*r/(dble(16.0)*scl)
     func_covy=(dble(1.0)+r*r/dble(8.))*r*bessk1(r)-r*r*bessk0(r)
!     !
           end if
! !c#1       
       else
          func_covy=dble(1.0)
! !c#1
       end if
      
  IF(indec.eq.5)THEN
     IF(r.gt.1.0d-10)THEN
        func_covy=c0*r2**omega2
     ELSE
        func_covy=0.0d0
     ENDIF  
  ENDIF  
 
! write(*,*)'in covy (indec,x,y,FUNC_COVY)',indec,x,y,func_covy
 
      
ENDFUNCTION
            
      
