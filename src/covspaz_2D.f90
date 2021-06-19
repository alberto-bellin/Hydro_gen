!
 subroutine cov_spaz(n1,n2,lag_x,lag_y,field)

!
!
!*************************************************************************
!
!          this subroutine computes the spatial covariance with
!         reference to the generic realization.
!         the field is stored in the vector a of lenght n1*n2
!         where nx and ny are the # of points along x and y
!         directions.
!========================================================================
!                    dummy parameters description:
!
!        n1    => # of grid point in longitudinal direction
!        n2    => # of grid points in transverse direction
!        ngen   => n1*n2 (total number of grid points)
!        lag_x => max lag in longitudinal direction 
!        lag_y => max lag in transverse direction 
!        field => real*8 vector containing the generated field stored
!                 line by line
!      covarx  => real*8 vector containing the longitudinal spatial
!                 covariance function.
!      covary  => real*8 vector containing the transverse spatial
!                 covariance function.
!    
!
!======================================================================= 
!
!

 use vars_common_hydrogen, only: dx,dy

implicit none


integer                            :: n1,n2,ngen,lag_x,lag_y,lag,lagmx,lagmin,np
integer                            :: i,j,k,iniz,fin

real                               :: mean
!real                               :: covarx(*),covary(*)
real*8                             :: field(n1*n2)
real                               :: scalex,scaley
real, allocatable,    dimension(:) :: covarx,covary



! compute the mean

 ngen=n1*n2

 
 mean=0.0d0
 mean = sum(field(:))
 mean= mean/ngen
 field=field-mean

      
!
! transverse covariance function
!
 allocate(covarx(lag_x),covary(lag_y))
 

 do lag=1,lag_y
    
    covary(lag)=0.0d0
    np=n2*(n1-lag+1)
      do j=1,n2
         do i=1,n1-lag+1
            iniz=(i-1)*n1+j
            fin=iniz+lag-1
            covary(lag) = covary(lag)+field(iniz)*field(fin)
         end do
      end do
    
    covary(lag)=covary(lag)/np
 end do


!
!  longitudinal covariance function
!

 do lag=1,lag_x

    covarx(lag)=0.0d0
    np=n1*(n2-lag+1)
      do i=1,n1
         do j=1,n2-lag+1
            iniz=(j-1)*n1+i
            fin=(j+lag-2)*n1+i
            covarx(lag) = covarx(lag)+field(iniz)*field(fin)
         end do
      end do
    
    covarx(lag)=covarx(lag)/np
 end do

       field =field+mean
 
       
        
        open(99,file='spatial_statistics.txt',form='formatted')
        
        WRITE(*,*)'      MEAN          VARIANCE  '
        write(*,*)mean,COVARX(1)
        WRITE(99,*)'      MEAN          VARIANCE  '
        write(99,*)mean,COVARX(1)
        write(99,*)
        
        write(99,*)'   lag   covx      covy     '
        write(99,*)'------------------------------------'

        lagmx=max(lag_x,lag_y)
        lagmin=min(lag_x,lag_y)
        
  
        if(lagmx-lagmin.eq.0) then

            do i=1,lagmx
                 write(99,102)i,covarx(i),covary(i)
            end do

        else

            do i=1,lagmx
                if(i.le.lagmin) then
                    write(99,102)i,covarx(i),covary(i)
                end if

                if(i.gt.lag_x.and.i.le.lag_y) then
                    write(99,103)i,covary(i)
                end if

                if(i.le.lag_x.and.i.gt.lag_y) then
                    write(99,108)i,covarx(i)
                end if
            end do

        end if

        close(99)
        
102     format(2x,i4,2x,2(f8.6,2x))
103     format(2x,i4,10x,f8.6)
108     format(2x,i4,2x,f8.6)
    
   deallocate(covarx,covary)

 end subroutine






 
