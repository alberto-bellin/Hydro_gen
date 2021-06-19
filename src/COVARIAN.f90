!
      module mod_covar
      contains
      
      subroutine covar(filecov,idim,ilrf,nn1,jfin,indec,dx,dy, &
     sclx,scly,c0,omega2,cov,cov1,cov2)  

     
     
     
!      use vars_common_hydrogen

!
!----- This subroutine compute the covariance matrix
!      or the semivariogram matrix when power law semivariograms
!      are employed. 
!      It uses the function func_covy(x,y) which gives the
!      covariance function or the semivariogram for the lag (x,y)
!
!      For kriging systems defined by the semivariogram the kriging
!      coefficients are computed with reference to a pseudo covariance
!      function defined as: C(x_x,x_j)=A-gamma(x_i,x_j), where A
!      is a positive constant greater than the greatest semivariogram used
!      in the kriging system (Journell and Huijbregts, Mining Geostatistics
!      pg 306).
!           This subroutine computes the variogram and the main program
!      computes the pseudo covariance employed in the kriging system
!

      use mod_func_covy
      
      implicit none
      
    
      
      character*30, intent(in)  :: filecov
      
      integer, intent(in) ::  idim
      integer, intent(in) ::  ilrf
      integer, intent(in) ::  nn1
      integer, intent(in) ::  jfin
      integer, intent(in) ::  indec
      
      real*8,  intent(in) :: dx
      real*8,  intent(in) :: dy
      real*8,  intent(in) :: sclx
      real*8,  intent(in) :: scly
      real*8,  intent(in) :: c0
      real*8,  intent(in) :: omega2
      
      real*8, dimension(idim,jfin),intent(out) :: cov
      real*8, dimension(3,4,ilrf),intent(out)  :: cov1
      real*8, dimension(10,ilrf),intent(out)   :: cov2
      
      
      real*8, parameter  :: zero=0.0d0
      real*8, parameter  :: two=2.0d0
      integer, parameter :: ilevmx=4
     
     
      integer  :: i,j,kk
      real*8   :: dx1,dy1,x,y
      
      
    
!      dimension cov(idim,*),cov1(3,4,ilevmx),cov2(10,ilevmx)
!      character*30 filecov
!      external covy
      
      

      
      if(ilrf.gt.ilevmx) then
          write(*,*)'increases the parameter ilevmx in the files'
          write(*,*)'covariance.f and hydro_gen.fto the following'
          write(*,*)'value:',ilrf
          stop
      end if

      if(indec.eq.0) then
!
!---- the covariance function is read from a file ---
!
      open(88,file=filecov)
     
      do i=1,nn1+1
      read(88,*)(cov(i,j),j=1,jfin)
      end do  
      
      do kk=1,ilrf
!----- level # 1 ------
       do i=1,2
      read(88,*)(cov1(i,j,kk),j=1,2)
       end do
       read(88,*)(cov1(3,j,kk),j=1,4)
!
!----- level # 2 -------
!
         do i=1,10
         read(88,*)cov2(i,kk)
           end do 
      end do

       close(88,status='keep')

      
       else

      do i=1,nn1+1
        y=dble(i-1)*dy/scly
      do j=1,jfin
        x=dble(j-1)*dx/sclx
        cov(i,j)=func_covy(x,y,c0,omega2,sclx,indec)
      end do
      end do
      
 
        dx1=dx/sclx
        dy1=dy/scly  


        
        do kk=1,ilrf
        
!
!------ level # 1------
!
        do i=1,2
        y=dble(i-1)*dy1
        do j=1,2
        x=dble(j-1)*dx1
        cov1(i,j,kk)=func_covy(x,y,c0,omega2,sclx,indec)
        end do
        end do
        cov1(3,1,kk)=func_covy(-dx1/two,-dy1/two,c0,omega2,sclx,indec)
        cov1(3,2,kk)=func_covy(dx1/two,-dy1/two,c0,omega2,sclx,indec)
        cov1(3,3,kk)=func_covy(-dx1/two,dy1/two,c0,omega2,sclx,indec)
        cov1(3,4,kk)=func_covy(dx1/two,dy1/two,c0,omega2,sclx,indec)
!
!---- level # 2
!

      cov2(1,kk)=func_covy(zero,zero,c0,omega2,sclx,indec)
      cov2(2,kk)=func_covy(dx1/two,-dy1/two,c0,omega2,sclx,indec)
      cov2(3,kk)=func_covy(dx1,zero,c0,omega2,sclx,indec)
      cov2(4,kk)=func_covy(dx1/two,dy1/two,c0,omega2,sclx,indec)
      cov2(5,kk)=func_covy(zero,dy1,c0,omega2,sclx,indec)
      cov2(6,kk)=func_covy(-dx1/two,dy1/two,c0,omega2,sclx,indec)
      cov2(7,kk)=func_covy(-dx1/two,zero,c0,omega2,sclx,indec)
      cov2(8,kk)=func_covy(zero,-dy1/two,c0,omega2,sclx,indec)
      cov2(9,kk)=func_covy(dx1/two,zero,c0,omega2,sclx,indec)
      cov2(10,kk)=func_covy(zero,dy1/two,c0,omega2,sclx,indec)
                
       dx1=dx1/two
       dy1=dy1/two
       end do
       

       end if

!      write(*,*)'------- covariance matrix -----'
!      do i=1,nn1+1
!      write(*,888)(cov(i,j),j=1,jfin)
!      end do
!         write(*,*)'level # 1'
!         do i=1,2
!       write(*,*)(cov1(i,j),j=1,2)
!         end do
!       write(*,*)(cov1(3,j),j=1,4)
!
!        write(*,*)'level # 2'
!         do i=1,10
!          write(*,*)cov2(i)
!         end do 
!888   format(10f8.5)
      return
      end subroutine
      end module
 
