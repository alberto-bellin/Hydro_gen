! program  coefcy.f
!
!
      module mod_krig
      contains
      
      subroutine krig(ilevel,nm,jfin,nrow,ncol,cov,vcx)
  
!
!----- ilevel: 0 -> coarse grid; 1-> first level of the refinement 2-> second level of the refinement
!      nm -> dimension of the kriging matrix
!      jfin -> number of nodes along the direction x (j) in the searching area
!      nrow -> number of rows of the covariance matrix
!      ncol -> number of columns of the covariance matrix
!      cov  -> covariance matrix
!      conditional variance
!

use vars_common_hydrogen,only:u1,itype

!      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      implicit none

      integer,intent(in)  :: jfin
      integer,intent(in)  :: nm
      integer,intent(in)  :: ilevel
      integer,intent(in)  :: nrow
      integer,intent(in)  :: ncol
      
      real*8,dimension(nrow,ncol) :: cov
      
      
      real*8,intent(out)  :: vcx
 
!      include 'dimension.blk' 
!      parameter(icond=1670,iwork=icond*(icond+1)/2)
 
 integer  :: i,j,k,i1,i2,n1,n2,ilag,jlag,ntot
 integer  :: info
 integer  :: iwork
 real*8   :: zero,one

 integer,allocatable,dimension(:)      ::iwk
 real*8,allocatable,dimension(:)       :: ap
 real*8,allocatable,dimension(:)       :: gama

 
 

!      DIMENSION ap(iwork),GAMA(icond),U1(*),cov(idm,*)
!        integer iwk(iwork)

        data zero,one/0.0d0,1.0d0/

! iwork=ijmax*(ijmax+1)/2
  iwork=nm*(nm+1)/2
  
 allocate(ap(iwork))
 allocate(iwk(iwork))
 allocate(gama(nm))
 
 

!
!
!-- compute the kriging matrix ----
!

	N1=nm-1
        k=0
	DO  j=1,N1
            i2=int((j-1)/jfin)                                  ! first node     
              do i=1,j
                 i1=int((I-1)/jfin)                             ! second node
                 ilag=i2-i1+1
                 jlag=iabs(j-i-(ilag-1)*jfin)+1
                 k=k+1
                 ap(k)=cov(ilag,jlag)                           ! only the upper triangle of the matrix is stored by columns
              end do
        end do


       
!
!
!---- compute the known term vector 
!
      if(ilevel.eq.0) then

          i2=int((nm-1)/jfin)
          DO  I=1,N1
          i1=int((i-1)/jfin)
          ilag=i2-i1+1
          jlag=iabs(nm-i-(ilag-1)*jfin)+1
          GAMA(i)=cov(ilag,jlag)
          u1(i)=gama(i)
          end do
       else

          do i=1,n1
          gama(i)=cov(3,i)
          u1(i)=gama(i)
          end do
       end if


      
       ntot=n1

       if(itype.eq.5) then
       do i=1,n1
          k=k+1
          ap(k)=one
       end do
       ap(k+1)=zero
       gama(nm)=one
       u1(nm)=gama(nm)
       ntot=nm
       end if
!
!--------------------------------------------------------------------
!         factors a double precision symmetric matrix stored in
!         packed form by elimination with symmetric pivoting.
!----------------------------------------------------------------------
!
      info=0
      call dspfa(ap,ntot,iwk,info)

       if(info.ne.0) then
           write(*,*)'error in matrix factorization'
           stop
        end if
!
!---------------------------------------------------------------------
!     solves the double precision symmetric system
!     a * x = b
!     using the factors computed by dspfa.
!---------------------------------------------------------------------

      call dspsl(ap,ntot,iwk,u1)

      VCX=dble(0.)


      DO I=1,ntot
      VCX=VCX+U1(i)*GAMA(i)
        end do
!       write(*,*)nd,vcx
!
!
 deallocate(ap)
 deallocate(gama)
 deallocate(iwk)
 
	RETURN
	END subroutine
	end module
