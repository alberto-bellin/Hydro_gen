program hydrogen

use vars_common_hydrogen,only:nseed,file1,iNP,np

!=============================================================================
!
!  H   H Y        Y DDD     RRRRR      OOOO      GGGGGGG  EEEEEE N         N
!  H   H  Y      Y  D  D    R    R    O    O     G        E      N N       N
!  H   H   Y    Y   D   D   R     R  O      O    G        E      N  N      N
!  H   H    Y  Y    D    D  R     R  O      O    G        E      N   N     N
!  HHHHH     Y      D    D  R    R   O      O    G        EEEE   N    N    N
!  H   H     Y      D    D  R  R     O      O    G   GGG  E      N     N   N
!  H   H     Y      D   D   R   R    O      O    G     G  E      N      N  N
!  H   H     Y      D  D    R    R    O    O     G     G  E      H       N N
!  H   H     Y      DDD     R     R    OOO  ____ GGGGGGG  EEEEEE H         N
!
!============================================================================
!
! Generation of normally ditributed random fields with a given covariance
! function
!                               BY:
!
!               ALBERTO BELLIN^1 AND YORAM RUBIN^2
!
!       1:   Dipartimento di Ingegneria Civile, Ambientale e Meccanica
!            Universita' di Trento
!            via Mesiano, 77, I-38050 Trento, Italy
!            phone:  +39 461 882620
!            fax:    +39 461 882672
!            e-mail: alberto@unitn.it
!
!       2:   Department of Civil Engineering
!            University of California, Berkeley
!            Berkeley, CA 94720, USA
!            phone:  +1 510 642 2282
!            fax:    +1 510 642 7476
!            e-mail: rubin@ce.berkeley.edu
!
!=============================================================================
!
! Software developed by Alberto Bellin  
! Universita' di Trento,
! Dipartimento di Ingegneria Civile, Ambientale e Meccanica
! 38050-I Mesiano di Povo, TRENTO
!
! e-mail:   Alberto.Bellin@unitn.it 
!
! Copyright: Alberto Bellin and Yoram Rubin.
!
!
!
! Summary:  This is the main program of the package aimed at the generation
!           of two dimensional random fields with assigned covariance
!           function.
!
!
!
! Package Version: 2.0 April 1997
!
!
! please cite the following paper in papers or reports that use
! the present software:
!
! Bellin A., Y. Rubin, Hydro_gen: A new random field generator
! for correlated properties, Stochastic Hydrology and Hydraulics,
! 10(4), 1996.
!
!
! Permission is  granted to anyone to use and modify this packages provided
! that:
! i) the authors are acknowledged by citing the abofe referenced paper;
! ii) the use in any kind of research or job will be cited in the relative
! papers or reports;
! iii) the use of the package is under the user responsability
! NO WARRANTY is given concerning bugs and errors.
! iv) The use or distribution must be free of charge.                 
!  v) The package uses the following libraries:
!     a) LINPACK by J. J. Dongarra, J. R. Bunch, C. B. Moler 
!        e G.W. Stewart, for the linear system solution
!     b) BLAS, for linear algebra
!     d) RANLIB by Barry W. Brown and James Lovato,
!        Department of Biomathematics, Box 237
!        the University of Texas, M.D. Anderson Cancer Center
!        1515 Holcombe Boulevard, Huston, TX 77030, for the generation
!        of independent normally distributed random numbers.
!     e) Numerical Recipes by W. H. Press, B. P. Flannery, S. A.
!        Teukolsky, W. T. Vetterling, for the function computing
!        the Bessel Function
!  
! Copyright conditions of the above referenced libraries are 
! extended to hydro_gen  
!
! Bug reports and hints are welcomed to the following e-mail address:
! alberto.bellin@unitn.it
!----------------------------------------------------------------------------
!              parameter dimensions description:
!
!  igrid1 => max number of grid points in y direction
!  igrid2 => max number of grid points in x direction
!            BOTH ARE FOR THE COARSE GRID
!  icond  => max number of kriging points
!  iptmx  => max dimension for the vectors storing the kriging
!            coefficients
!            the coefficients are stored in the following order:
!            side #1: kriging area exiting form the field border
!                      x=0
!            side #2: kriging area exiting from the field border
!                     x=lx
!            side #3: kriging area exiting from the field border
!                     y=0
!            corners: kriging area exiting from the two couple
!                     of borders y=0 ,x=0 or y=0, x=lx.
!
!  iptmxcv=> max dimension for the  vectors containing the
!            conditional variances relatives to points referenced
!            in the previous entry.
!            the coefficients are stored in the following order:
!            side #1: kriging area exiting form the field border
!                      x=0
!            side #2: kriging area exiting from the field border
!                     x=lx
!            side #3: kriging area exiting from the field border
!                     y=0
!            corners: kriging area exiting from the two couple
!                     of border y=0 ,x=0 or y=0, x=lx.
!
!  idim   => max dimension of the vector storing the covariance
!            function to be reproduced
!
!
!
!
!
!---------------------------------------------------------------------
!       cond(i,j) -> generated field value at position i,j
!
!
!       u1sid(I) -> vector of the kriging coefficients
!                   which are used in the Monte Carlo generations
!
!       vcxsid(i) -> vector of  the conditional variances
!cdx,dy
!
!--------------------------------------------------------------------------


 !read the input parameters
 implicit none


 !REAL:: T1,T2

 !CALL CPU_TIME(T1)
 
 call hydrogen2D_input
 write(*,*) '---------------------------'
 write(*,*) 'hydrogen2D_input has ended '
 write(*,*) '---------------------------'

 open(18,file=file1,status='unknown')
 open(20,file='stats.out',status='unknown')
 
 call coeff_2D
 
 write(*,*) '---------------------------'
 write(*,*) 'coeff_2D has ended '
 write(*,*) '---------------------------'
 
 !variabile file1 riga 197 hydrogen_input.90
 !read(88)c,beta  riga 184 variabile cn non dichiarata. esiste?
 !
 !
 !covartype vettore di dimensione 5, vedi mod_commonvars
 !
 !   

!----------------START MONTE CARLO SIMULATION--------------
!
!  NP  Number of Monte Carlo Simulations.
!  i    row indicator
!  j    column indicator
!

      WRITE(*,*)
      WRITE(*,*)' kriging coefficients computed; starting Monte Carlo..'

        
!  iprint=0
MC_loop: DO  iNP=1,NP
 

 !       initialize the random number sequence
!       original code in f77
!       call setall(nseed,nseed1) 


!       new call for ranlib in f90
!        call set_initial_seed (nseed, nseed1 )               
        
!       
!        write(*,*)'SEEDS:',nseed,nseed1
!       write(*,*)
 
 call gen_field2D
 
 
  end do MC_loop

! CALL CPU_TIME(T2)
! WRITE(*,*) 'Elapsed time = ', T2-T1
end program
