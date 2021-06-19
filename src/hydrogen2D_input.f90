subroutine hydrogen2D_input
 use vars_common_hydrogen
 implicit none
 real*8        :: xspsug,xspasug,yspsug
 character*30  :: filecoef
      
      OPEN(unit=1984, file='hydrogen_input.txt',action='read')
      read(1984,*)nseed
      read(1984,*)dx,dy
      read(1984,*)np
      read(1984,*)lx,ly
      read(1984,*)imark
      read(1984,*)itype
      read(1984,*)sigy
      read(1984,*)cond10
      read(1984,*)sclx,scly
      read(1984,*)xsp,ysp
      read(1984,*)xspa
      read(1984,10)filecoef  
      read(1984,*)ilevref
      read(1984,10) file1
      read(1984,*)iformat
      read(1984,*)flag_stat
      CLOSE(1984)
      write(*,*)'flag_stat:',flag_stat
      
      open(88,file=filecoef,form='unformatted',status='unknown')      
!      write(*,*) 'end of reading'
      
      covartype(1) = 'Exponential '
      covartype(2) = 'Gaussian '
      covartype(3) = 'Whittle '
      
!      ,'Gaussian ','Whittle ','da inserire','da inserire'/)
      
!! SEB      
!!      write(*,*)'Enter the seed number as an integer'
!!      write(*,*)'between 2 and 2147483398:'
!!      read(*,*)nseed
!!!               nseed1 = nseed-1
!!
!!
!! 
!!
!!      
!! 
!!
!!!-------   INPUT ------------------
!!!      DX: grid dimension along x
!!!      DY: grid dimension along y
!!!      lx & ly: field dimensions along x and y respectively
!!!      sigy: log-conductivity variance
!!!      COND10: mean log-conductivity
!!!      xsp,ysp: dimensions of the area used for the conditioning
!!!
!!!-------------------------------------------------------------
!!        write(*,*)'Enter final grid spacing dx,dy'
!!        read(*,*)dx,dy
!!        write(*,'(2x,A,2F6.2)')'dx and dy',dx,dy
!!        write(*,*)'Enter number of Monte Carlo realizations'
!!        read(*,*)np
!!        write(*,*)'# of monte carlo iterations:',np
!!        write(*,*)'Enter field dimensions'
!!        read(*,*)lx,ly
!!        write(*,'(2X,A,2F8.2)')'field dimensions:',lx,ly
!!        write(*,*)'To read the interpolation "kriging" coefficients'
!!        write(*,*)'from a file enter [1]; if not,'
!!        write(*,*)'enter any other integer:'
!!        read(*,*)imark
!!        write(*,*)'Enter the covariance type:'
!!        write(*,*)'itype=0 ==> discrete covariance function (from file)'
!!        write(*,*)'itype=1 ==> exponential'
!!        write(*,*)'itype=2 ==> Gaussian'
!!        write(*,*)'itype=3 ==> Whittle'
!!        write(*,*)'itype=4 ==> Mizell (B)'
!!        write(*,*)'itype=5 ==> power law semivariogram'
!!        write(*,*)'(self similar field)'
!!        read(*,*)itype
!!
!!        write(*,101)covartype(itype)
!!101     format('Covariance type:',a11)
!!
!!        write(*,*)'covariance of type:',itype 
!!        
!!        if(itype.ne.5) then          
!!           write(*,'(2X,A,F8.2)')'Enter the variance'
!!           read(*,*)sigy
!!           write(*,*)'variance:',sigy
!!        end if
!!      
!!        if(itype.eq.5) then
!!          write(*,*)'do you want to force the mean to be constant in'
!!          write(*,*)'each realization? [Y/N]'
!!          read(*,'(A1)')answ
!!          if(answ.eq.'Y') then
!!             write(*,*)'enter the mean:'
!!             read(*,*)fixcond
!!             cond10=0.0d0
!!             write(*,'(2X,A,F8.2)')'mean:',fixcond
!!          end if                                 
!!        else
!!          write(*,*)'Enter the mean '
!!          read(*,*)cond10
!!          write(*,'(2X,A,F8.2)')'mean:',cond10        
!!        end if
!!        
!!        
!!        if(imark.ne.1)then !if_imark
!!
!!        if(itype.eq.1) then
!!         write(*,*)'Enter integral scales in x and y directions:'
!!         read(*,*)sclx,scly
!!         write(*,'(2X,A,2F8.2)')'integral scales:',sclx,scly
!!         xspsug = DBLE(3.0)*sclx
!!         yspsug = DBLE(4.0)*scly
!!         xspasug = sclx
!!        end if
!!
!!        if(itype.eq.2) then
!!         write(*,*)'Enter correlation lengths in x and y directions:'
!!         read(*,*)sclx,scly
!!         write(*,'(2X,A,2F8.2)')'correlation lengths:',sclx,scly
!!         xspsug = DBLE(5.0)*sclx
!!         yspsug = DBLE(6.0)*scly
!!         xspasug = sclx
!!        end if
!!
!!      if(itype.eq.3.or.itype.eq.4) then
!!         write(*,*)'Enter the integral scale (isotropic)'
!!         read(*,*)sclx
!!         write(*,'(2X,A,F8.2)')'integral scale (isotropic):',sclx
!!         scly=sclx
!!        if (itype .eq. 3) then
!!           xspsug = DBLE(3.0)*sclx
!!           yspsug = DBLE(4.0)*scly
!!           xspasug = sclx
!!        else
!!           xspsug = DBLE(5.0)*sclx
!!           yspsug = DBLE(6.0)*scly
!!           xspasug = sclx
!!        end if
!!      end if
!!
!!        if(itype.eq.5) then
!!          write(*,*)'Enter the reference scales in x and y directions:'
!!          read(*,*)sclx,scly
!!          write(*,'(2X,A,2F8.2)')'reference  scales:',xspasug,sclx,scly
!!          write(*,*)'enter the cofficients a and beta of the power'
!!          write(*,*)'semivariogram ( g=a r^beta)'
!!          read(*,*) beta
!!        end if                                                
!!                                             
!!        if(itype.ne.5) then
!!                                       
!!        write(*,*)'search neighborhood dimensions:'
!!        write(*,*)'--------           -'
!!        write(*,*)'|      |           |'
!!        write(*,*)'|      --------   ysp '
!!        write(*,*)'|             |    |'
!!        write(*,*)'|             |    |'
!!        write(*,*)'---------------    -'
!!        write(*,*)'|--xsp-|-xspa-|'
!!        write(*,*)
!!        write(*,*)'Suggested values:' 
!!        write(*,'(2X,A,F8.3)')'xsp = ',xspsug
!!        write(*,'(2X,A,F8.3)')'ysp =',yspsug
!!        write(*,'(2X,A,F8.3)')'xspa =',xspasug
!!        write(*,*)
!!
!!        if (xspsug+xspasug .GT. lx-3*sclx) then
!!           write(*,*)'WARNING:'
!!           write(*,*)'Suggested xsp/xspa is too large.' 
!!           write(*,*)'Enter xsp and xspa such that their sum'
!!           write(*,*)'is smaller than:',lx-3*sclx
!!           write(*,*)'or, stop program and begin again with an'
!!           write(*,*)'x dimension larger than:',xspsug+xspasug+3*sclx
!!        end if
!!
!!        if (yspsug .GT. ly-3*scly) then
!!           write(*,*)'WARNING:'
!!           write(*,*)'Suggested ysp is too large.' 
!!           write(*,*)'Enter an ysp smaller than:',ly-3*scly
!!           write(*,*)'or, stop program and begin again with a'
!!           write(*,*)'y dimension larger than:',yspsug+3*scly
!!        end if
!!        write(*,*)
!!
!!        write(*,*)'Enter xsp, ysp'
!!        read(*,*)xsp,ysp
!!        write(*,'(2X,A,2F8.3)')'xsp and ysp:',xsp,ysp
!!        write(*,*)'Enter xspa'
!!        read(*,*)xspa
!!        write(*,'(2X,A,F8.3)')'xspa:',xspa
!!        end if
!!
!!      end if!end_if_imark
!!      
!!      write(*,*)'Enter the file name for the file storing the'
!!      write(*,*)'interpolation coefficients:'
!!      read(*,10)filecoef  
!!      open(88,file=filecoef,form='unformatted',status='unknown')
!!
!!      if(itype.le.4) then
!!        if(imark.eq.1) then
!!           read(88)xsp,ysp
!!           read(88)xspa 
!!        else
!!           write(88)xsp,ysp
!!           write(88)xspa
!!        end if        
!!      end if
!!      
!!     if(itype.eq.5) then
!!         if(imark.eq.1) then
!!            read(88)c,beta
!!         else
!!            write(88)c,beta
!!         end if
!!      end if
!!      !         
!!      write(*,*)'enter the number of refinement levels:'
!!      write(*,*)'[0] ==> No refinement'
!!      write(*,*)'[n] ==> refinement at n levels'
!!      read(*,*)ilevref
!!      
!!      if(ilevref.gt.4) then
!!          write(*,*)'Number of refinements set to 4'
!!          ilevref=4
!!      end if
!!      
!!      !
!!      write(*,*)'Enter the name  of the output file for the'
!!      write(*,*)'replicates [max 30 characters]:'
!!      read(*,10) file1
!!      write(*,*)'Enter the format of the output file:'
!!      write(*,*)'type 1 for 3 columns x y z'
!!      write(*,*)'type 0 for the matrix format (only z values'
!!      write(*,*)'on a regular grid)'
!!      read(*,*)iformat
!!      !
!! SEB
10    format(a30)

end subroutine

