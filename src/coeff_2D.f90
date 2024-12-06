subroutine coeff_2D
use vars_common_hydrogen
use mod_krig
use mod_covar
use mod_coefl2

implicit none
integer :: kin,nx,ny
integer :: jfin,jjfin,jfin1,nlgs,ierror,jend,jstart
integer :: iiniz,jiniz,inum
integer :: idec
integer :: idim

integer :: i,j,icol,kk,iline,k,ii,jj
integer :: ijcountm1,ijcount,kcount,jjfin1
integer :: igrid1,igrid2,icond

real*8  :: betad2,fd
real*8  :: spmx,xlgs,semivar
real*8  :: betaxlgs
real*8,parameter  :: one=1.0d0


real*8,allocatable,dimension(:)     :: r
real*8,allocatable,dimension(:,:)   :: cond 

 
!
! ---------------------------------------------------------
!     SD1: S.D. of the conductivity
!----------------------------------------------------------
        SD1=DSQRT(SIGY)
!
!2000    format('particle#: ',i5)
!
!-----------------------------------------------------------
!       N1: # of grid points along ynlgs
!       N2: # of grid points along x
!       NN1: # of previous points along y used for the conditioning
!       NN2: # of previous points along x used for the conditioning
!       NN2A: # of following points along x used for the conditioning
!                       


      if(ilevref.eq.0) then
!          iref=2
          idiv=1
      else
!          iref=1
          idiv=2**ilevref
      end if

          ddx=dx*dble(idiv)
          ddy=dy*dble(idiv)
          
      n(1)=int((ly+dy/10.)/dy)+1
      n(2)=int((lx+dx/10.)/dx)+1  
       
      kin=idiv  
      n1in=1
      do kk=1,ilevref
         kin=kin/2
         n1in=n1in+kin
      end do  
      n2in=n1in
                  
      n1fin=n1in+n(1)-1
      n2fin=n2in+n(2)-1
      
      nx=n2fin-n2in+1
      ny=n1fin-n1in+1   
      nxy=nx*ny
            
      if(ilevref.gt.0) then
         n1tot=n1fin+idiv+mod(n(1),2)
         n2tot=n2fin+idiv+mod(n(2),2)
         n1=int(n1tot/idiv)+1
         n2=int(n2tot/idiv)+1
      else
         n1tot=n(1)
         n2tot=n(2)
         n1=n1tot
         n2=n2tot
      end if
      

!      ngen=nx*ny
!      if(ilevref.gt.0) then
!      ngen=ngen+ilevref*(nx+ny+2*ilevref-2)+(n1+n2-2)*2  
!      end if
       ngen=((n1-1)*idiv+1)*((n2-1)*idiv+1)
      
      
!      write(*,*)'ngen:',ngen


!        igrid=2*igrid1
!        N2=INT(LX/DX)+1
!        N1=INT(LY/DY)+1
!
!        n2d2=int((2*n2-1)/2)
!        n1d2=int((2*n1-1)/2)
!
!            ilast=2*n1-1
!            jlast=2*n2-1
!        if(iref.eq.2) then
!         isrt=0
!          else
!         isrt=1
!          end if

!        if(iref.eq.1) then
!        n1n2=(2*n1-1)*(2*n2-1)
!        n1n2=n1n2-4*(n1+n2-2)
!        else
!        n1n2=n1*n2
!        end if
!        rn1n2=dble(n1n2)

!        NN2=int(xsp/dx)
!        NN1=int(ysp/dy)xlgs
!        NN2A=int(xspa/dx)
!        jfin=nn2+nn2a+1spmx
! 
!
!------- grid spacing for the coarse grid:
!

 

        if(itype.ne.5) then
          nn2=int(xsp/ddx)
          nn2a=int(xspa/ddx)
          nn1=int(ysp/ddy)
        else
!
!---- Power law semivariogram------
!          
          NN2=n2-1
          NN1=n1-1
          NN2A=0
        end if 

        jfin=nn2+nn2a+1
        idim=nn1+1
                  
      if(itype.eq.5) then   
          betad2=beta/dble(2.)
          fd=dble(3)-betad2
       
          write(*,*)'the semivariogram is computed with reference to'
          write(*,*)'the following parameters:'
          write(*,*)'C=',c
          write(*,*)'beta=',betaxlgs
          write(*,*)'FRACTAL DIMENSION:',fd
          write(*,*)'the employed semivariogram is stored in the file:'
          write(*,*)'semivariog.out'
          open(8,file='semivariog.out')
          spmx=dsqrt(lx*lx+ly*ly)
          nlgs=int(spmx/ddx)+1
          do i=1,nlgs
          xlgs=(i-1)*ddx
          semivar=c*xlgs**beta
          write(8,*)xlgs,semivar
          end do
          close(8,status='keep')
       end if

!! SEB serve inizializzazione altrimenti si attiva lo stop nell'if che segue                         
        ierror=0
!! SEB
!        if(n1tot.gt.igrid1) then
!            write(*,*)'DIMENSIONING ERROR: the parameter igrid1'
!            write(*,*)'in the file dimension.blk is not properly'
!            write(*,*)'assigned. Suggested value:',n1tot+1
!            ierror=ierror+1
!         end if
! 
!         if(n2tot.gt.igrid2) then
!            write(*,*)'DIMENSIONING ERROR: the parameter igrid2'
!            write(*,*)'in the file dimesnion.blk is not properly'
!            write(*,*)'assigned. Suggested value:',n2tot+1
!            ierror=ierror+1
!         end if
 
        ijcount=(nn2+nn2a+1)*nn1+nn2

!         if(ijcount.gt.icond) then
!            write(*,*)'DIMENSIONING ERROR: the parameter icond'
!            write(*,*)'in the file dimension.blk is not properly'
!            write(*,*)'assigned. Suggested value:',ijcount
!            ierror=ierror+1
!         end if

!! SEB decommentata
        if(ierror.gt.0) stop

         write(18,101)covartype(itype)
         write(18,102)np
102     format('This file stores',i4,'independent replicate(s)')
        
        if(iformat.eq.0) then
          write(18,*)'the data are stored in the matrix format with:'
          write(18,67)nx,ny        
        else
          write(18,*)'the data are stored in column format (x y z)'
          write(*,*)'total number of data:',nx*ny

        end if


67     format(i3,1x,'lines and',1x,i3,1x,'columns')
      write(18,*)'grid size:',dx,dy
      write(18,*)'---------------------------------------------------'  
          
      write(20,*)'---------------------------------------------------'
      write(20,*)'replicate N.           mean        variance'
      write(20,*)'-------------------------------------------'
         
!!SEB
!!SEB aggiunto formato 101
101     format('Covariance type:',a11)
!! SEB
!
!------ compute the covariance matrix --
!
!      itype=0   discrete covariance function (it is read from a file)
!      itype=1   exponential covariance function (defined by a function)

!
       if(itype.eq.0) then
          write(*,*)'the covariance function is read from the file:',filecov
          read(*,442)filecov
       end if
       442 format(a30)


!-----------------------------------------------
!


!
!CCCCCCCC    RUN TO IDENTIFY THE DIMENSION OF VECTORS STORIG THE KRIGING COEFFICIENTS AND THE CONDITIONAL VARIANCES
!
!

	if(itype.eq.5) then

        write(*,'(2x,A)')'line   nvectpos   kkcv'
      jfin=n2 
      i=1
      do j=2,n2      
              ijcountm1=j-1
             do kcount=1,ijcountm1
                nvectpos=nvectpos+1
             end do             
 
             kkcv=kkcv+1    
      end do 
         write(*,'(2x,i3,i8,i8)')i,nvectpos,kkcv  

      do i=2,n1            
         do j=1,n2  
             ijcount=(i-1)*n2+j  
             ijcountm1=ijcount-1
!           write(*,*)i,j
             
             do kcount=1,ijcountm1
                nvectpos=nvectpos+1
             end do    
             kkcv=kkcv+1
             
         end do 
!         write(*,*)'line:',i,'done!'
         write(*,'(2x,i3,i8,i8)')i,nvectpos,kkcv
      end do   


	else

!
!-------- side #1
!
      kkcv=0
      nvectpos=0
      do icol=1,nn2
      ijcount=(nn2a+icol)*nn1+icol

      kkcv=kkcv+1
      jfin=nn2a+icol
      
      do kcount=1,ijcount-1
         nvectpos=nvectpos+1
      end do

      
      end do

!
!-------- side #2
!
      do icol=1,nn2a
      jfin=nn2+nn2a+1-icol
      ijcount=jfin*nn1+nn2+1
      kkcv=kkcv+1

      do kcount=1,ijcount-1
         nvectpos=nvectpos+1
      end do               
 
      end do      

!
!------- side #3
!
      jjfin=nn2+nn2a+1
      jfin1=nn2+1
                do iline=1,nn1
      if(iline.eq.1)  then
             jfin=jfin1
          else
             jfin=jjfin
          end if

        ijcount=jfin*(iline-1)+jfin1

       kkcv=kkcv+1      
      do kcount=1,ijcount-1
         nvectpos=nvectpos+1
      end do  
      end do


!
!----- corners
!

        DO 1110 i=1,nn1
        jstart=1
        jend=nn2

1120     DO 1100 J=jstart,jend
!
!
!
!----- select the points used for the kriging--------
       IINIZ=I-NN1
       JINIZ=J-NN2
       idec=1
        jfin=j+nn2a
         if(jiniz.ge.1.and.jfin.le.n2)  goto 1100
       iiniz=max(iiniz,1)
       jiniz=max(jiniz,1)
       jfin=min(jfin,n2)
        jfin1=j
        if(i.eq.1) jfin=j
!
!----- select the points used for the kriging--------
!
        jjfin=jfin-jiniz+1
        jjfin1=jfin1-jiniz+1
        inum=I-IINIZ
        ijcount=jjfin*inum+jjfin1

                if(ijcount.eq.1)  go to 1100

!
!------------- compute the kriging coefficients -------
      kkcv=kkcv+1
 
      do kcount=1,ijcount-1
         nvectpos=nvectpos+1
      end do
      

1100      continue
      if(jstart.eq.1) then
      jstart=n2-nn2a+1
      jend=n2
        go to 1120
          end if
1110   continue

!
!----------- KRIGING COEFFICIENTS FOR THE LARGER CONDITIONING AREA
!
      jfin=nn2+nn2a+1
      ijcount=nn1*jfin+nn2+1
      ijmax=ijcount

!------- compute the kriging coefficients-------
!      write(*,*)'compute the kriging coefficients for the larger area'
      kkcv=kkcv+1
     
      do k=1,ijcount-1

         nvectpos=nvectpos+1
      end do   
      iptcvst=kkcv
      

      end if       


!       end if
!        write(*,*)'total # of vector positions:',nvectpos,kkcv



 !
 !CCCCCCC    END OF THE RUN FOR DIMENSIONING
 !
 allocate (u1sid(nvectpos))
 allocate (vcxsid(kkcv))
 allocate (u1l1(4,ilevref))
 allocate (u1l2(4,ilevref)) 
 allocate (vcxl1(ilevref))
 allocate (vcxl2(ilevref))
 allocate (r(ngen))
 allocate (cond(n1tot+1,n2tot+1))
 allocate (u1(ijmax))
 allocate (cov1(3,4,ilevref))
 allocate (cov2(10,ilevref))
 allocate (cov(idim,jfin))
 allocate (cov1a(3,4))                  !! SEB modificato da 2,4 a 3,4
 allocate (cov2a(10))

nvectpos=0
kkcv=0
ierror = 0

if(imark.ne.1) then
! 
!       
!          ierror=max(nn1+1,jfin)
!          if(ierror.gt.idim) then
!          write(*,*)'DIMENSIONING ERROR: the parameter idim'
!          write(*,*)'in the file hydro_gen.f is not properly'
!          write(*,*)'assigned. Suggested value:',ierror
!          stop
!       end if
!
     
      call covar(filecov,idim,ilevref,nn1,jfin,itype,ddx,ddy,sclx,scly,c,betad2,cov,cov1,cov2)
      
end if
 

      
      if(itype.eq.5) then
!
!----  self-similar RF---------
!
      if(imark.ne.1) then
!
!---- compute the kriging coefficients
!                               
      write(*,'(2x,A)')'line   nvectpos   kkcv'
      jfin=n2 
      i=1
      do j=2,n2      
              ijcountm1=j-1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           call krig(0,j,idim,j,itype,cov,u1,vcx) 
           
              call krig(0,j,j,idim,jfin,cov,vcx) 
              write(88)(u1(kcount),kcount=1,ijcountm1)
              write(88)dsqrt(vcx)
             do kcount=1,ijcountm1
                nvectpos=nvectpos+1
                u1sid(nvectpos)=u1(kcount)
             end do             
 
             kkcv=kkcv+1
             vcxsid(kkcv)=dsqrt(vcx)       
      end do 
         write(*,'(2x,i3,i8,i8)')i,nvectpos,kkcv  

      do i=2,n1            
         do j=1,n2  
             ijcount=(i-1)*n2+j  
             ijcountm1=ijcount-1
!            write(*,*)i,j

             call krig(0,ijcount,n2,idim,jfin,cov,vcx) 
             write(88)(u1(kcount),kcount=1,ijcountm1)
             write(88)dsqrt(vcx)
             
             write(*,*)(u1(kcount),kcount=1,ijcountm1)
             write(*,*)dsqrt(vcx)
             
             do kcount=1,ijcountm1
                nvectpos=nvectpos+1
                u1sid(nvectpos)=u1(kcount)
             end do    
             kkcv=kkcv+1
             vcxsid(kkcv)=dsqrt(vcx)

         end do 
!         write(*,*)'line:',i,'done!'
         write(*,'(2x,i3,i8,i8)')i,nvectpos,kkcv
      end do   

      else
!
!----- read the kriging coefficients
!
       i=1 
      write(*,'(2x,A)')'line   nvectpos   kkcv'       
       do j=2,n2
          kkcv=kkcv+1
          read(88)(u1sid(k+nvectpos),k=1,j-1)
          read(88)vcxsid(kkcv)
          nvectpos=nvectpos+j-1
       end do

      write(*,'(2x,i3,i8,i8)')i,nvectpos,kkcv     
       do i=2,n1
          do j=1,n2
             ijcount=(i-1)*n2+j-1
             kkcv=kkcv+1
             read(88)(u1sid(k+nvectpos),k=1,ijcount)
             read(88)vcxsid(kkcv)
             nvectpos=nvectpos+ijcount
          end do
!        write(*,*)'line:',i,'done!'                    
         write(*,'(2x,i3,i8,i8)')i,nvectpos,kkcv          
       end do

      end if           
      else
!
!
!----- compute the kriging coefficients for kriging area that exits
!      from the side x=0.0
!



       write(*,*)'side #1'


      do icol=1,nn2
      ijcount=(nn2a+icol)*nn1+icol

!---------- compute the kriging coefficients -------
!
      kkcv=kkcv+1
      if(imark.ne.1) then
      jfin=nn2a+icol
      
      call krig(0,ijcount,jfin,idim,jfin,cov,vcxsid(kkcv))
      write(88)(u1(j),j=1,ijcount-1)
      do kcount=1,ijcount-1
         nvectpos=nvectpos+1
         u1sid(nvectpos)=u1(kcount)
      end do
      vcxsid(kkcv)=dsqrt(one-vcxsid(kkcv))
      write(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv) 

      else
      read(88)(u1sid(nvectpos+kcount),kcount=1,ijcount-1)
      nvectpos=nvectpos+ijcount-1
      read(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      end if
      end do
!
!----- kriging coefficients for kriging area that exits
!      from the side x=lx
!
!------- iptkr2  pointer for the kriging coefficients of side #2
!                the vector position from iptkr2 to iptkr3-1 contain
!                the kriging coefficients of side #2.
!                The positions between 1 and iptkr2 contain
!                the coefficients of side #1
!        iptcv2  pointer for the conditional variances of side #2
!                the vector position from iptcv2 to iptcv3-1 contain
!                the conditional variances of side #2.
!                The positions between 1 and iptcv2 contain
!                the conditional variances of side #1
!      write(88)nvectpos,kkcv
       iptkr2=nvectpos
       iptcv2=kkcv

      write(*,*)'# of vector positions for the side #1:',iptkr2,kkcv

      write(*,*)'side #2'

      do icol=1,nn2a
      jfin=nn2+nn2a+1-icol
      ijcount=jfin*nn1+nn2+1
!------- compute the kriging coefficients--------------
      kkcv=kkcv+1
      if(imark.ne.1) then
      call krig(0,ijcount,jfin,idim,jfin,cov,vcxsid(kkcv))
      
      write(88)(u1(kcount),kcount=1,ijcount-1)
      vcxsid(kkcv)=dsqrt(one-vcxsid(kkcv))
      write(88)vcxsid(kkcv)      
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      do kcount=1,ijcount-1
         nvectpos=nvectpos+1
         u1sid(nvectpos)=u1(kcount)
      end do               
 
      else
      read(88)(u1sid(nvectpos+kcount),kcount=1,ijcount-1)
      nvectpos=nvectpos+ijcount-1
      read(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      end if
      end do
!
!------- iptkr3  pointer for the kriging coefficients of side #3
!                the vector position from iptkr3 to iptkrc-1 contain
!                the kriging coefficients of side #3.
!        iptcv3  pointer for the conditional variances of side #3
!                the vector position from iptcv3 to iptcvc-1 contain
!                the conditional variances of side #3.
!      write(88)nvectpos,kkcv     
      iptkr3=nvectpos
      iptcv3=kkcv

       write(*,*)'# of vector positions for the side #2:',nvectpos,kkcv

!
!------ kriging coefficients along the line y=0.0
!
        write(*,*)'side #3'

      jjfin=nn2+nn2a+1
      jfin1=nn2+1
                do iline=1,nn1
      if(iline.eq.1)  then
             jfin=jfin1
          else
             jfin=jjfin
          end if

        ijcount=jfin*(iline-1)+jfin1
!----- compute the kriging coefficients ------
       kkcv=kkcv+1
      if(imark.ne.1) then
      call krig(0,ijcount,jfin,idim,jfin,cov,vcxsid(kkcv))
      write(88)(u1(kcount),kcount=1,ijcount-1)
      vcxsid(kkcv)=dsqrt(one-vcxsid(kkcv))
      write(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      do kcount=1,ijcount-1
         nvectpos=nvectpos+1
         u1sid(nvectpos)=u1(kcount)
      end do  
   
      else
      read(88)(u1sid(nvectpos+kcount),kcount=1,ijcount-1)
      nvectpos=nvectpos+ijcount-1
      read(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      end if
      end do

!
!------- iptkrc  pointer for the kriging coefficients of corners
!                the vector position from iptkrc to the iptkrs contain
!                the kriging coefficients of corners.
!        iptcvc  pointer for the conditional variances of side #3
!                the vector position from iptcvc to the end contain
!                the conditional variances of corners.
       
!       write(88)nvectpos,kkcv   
       iptkrc=nvectpos
       iptcvc=kkcv

        write(*,*)'# of vector positions for the side #3:',nvectpos,kkcv
!
!----- compute the remaining kriging coefficients----
!
        write(*,*)'corners:'


        DO 111 I=1,NN1
        jstart=1
        jend=nn2

112     DO 110 J=jstart,jend
!
!
!
!----- select the points used for the kriging--------
       IINIZ=I-NN1
       JINIZ=J-NN2
       idec=1
        jfin=j+nn2a
         if(jiniz.ge.1.and.jfin.le.n2)  goto 110
       iiniz=max(iiniz,1)
       jiniz=max(jiniz,1)
       jfin=min(jfin,n2)
        jfin1=j
        if(i.eq.1) jfin=j
!
!----- select the points used for the kriging--------
!
        jjfin=jfin-jiniz+1
        jjfin1=jfin1-jiniz+1
        inum=I-IINIZ
        ijcount=jjfin*inum+jjfin1

                if(ijcount.eq.1)  go to 110

!
!------------- compute the kriging coefficients -------
      kkcv=kkcv+1
      if(imark.ne.1) then
      call krig(0,ijcount,jjfin,idim,jfin,cov,vcxsid(kkcv))
      write(88)(u1(kcount),kcount=1,ijcount-1)
      vcxsid(kkcv)=dsqrt(one-vcxsid(kkcv))
      write(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      do kcount=1,ijcount-1
         nvectpos=nvectpos+1
         u1sid(nvectpos)=u1(kcount)
      end do
      else
      read(88)(u1sid(nvectpos+kcount),kcount=1,ijcount-1)
      nvectpos=nvectpos+ijcount-1
      read(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      end if

110         continue
      if(jstart.eq.1) then
      jstart=n2-nn2a+1
      jend=n2
        go to 112
          end if
111   continue
!
!------- iptkrs  pointer for the kriging coefficients of the larger
!                kriging area
!                the vector position from iptkrs to the end of
!                vectors contain the kriging coefficients of corners.
!      write(88)nvectpos
      iptkrs=nvectpos
       write(*,*)'# of vector positions for the corners:',nvectpos,kkcv
!
!----------- KRIGING COEFFICIENTS FOR THE LARGER CONDITIONING AREA
!
      jfin=nn2+nn2a+1
      ijcount=nn1*jfin+nn2+1
      ijmax=ijcount

!------- compute the kriging coefficients-------
!      write(*,*)'compute the kriging coefficients for the larger area'
      kkcv=kkcv+1
      if(imark.ne.1) then
      call krig(0,ijcount,jfin,idim,jfin,cov,vcxsid(kkcv))
      write(88)(u1(j),j=1,ijcount-1)
      vcxsid(kkcv)=dsqrt(one-vcxsid(kkcv))
      write(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      do k=1,ijcount-1
         nvectpos=nvectpos+1
         u1sid(nvectpos)=u1(k)
      end do   
      iptcvst=kkcv
      else
      read(88)(u1sid(nvectpos+j),j=1,ijcount-1)
      nvectpos=nvectpos+ijcount-1
      read(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)   
      iptcvst=kkcv
      end if 

      end if       

!88       format(6i8)
       write(*,*)'total # of vector positions:',nvectpos,kkcv

!       ierror=0
!        if(iptmx.lt.nvectpos) then
!         write(*,*)'ERROR IN DIMENSIONING THE VECTOR STORING'
!         write(*,*)'THE INTERPOLATION COEFFICIENTS:'
!         write(*,*)
!         write(*,*)'ACTION: set the parameter iptmx to:',nvectpos
!         write(*,*)'in the file hydro_gen.f and recompile the'
!         write(*,*)'package. As an alternative, reduce the dimensions'
!         write(*,*)'of the larger search neighborhood or increase the'
!         write(*,*)'coarse grid spacing.'
!         ierror=ierror+1
!        end if

!        if(iptmxcv.lt.kkcv) then
!         write(*,*)'ERROR IN DIMENSIONING THE VECTOR STORING'
!         write(*,*)'THE CONDITIONAL VARIANCES:'
!         write(*,*)
!         write(*,*)'ACTION: set the parameter iptmxcv to:',kkcv
!         write(*,*)'in the file hydrogen.dim and compile again the'
!         write(*,*)'package. As an alternative reduce the dimensions'
!         write(*,*)'of the larger search neighborhood or increase the'
!         write(*,*)'coarse grid spacings'
!         ierror=ierror+1
!        end if
! 
!        if (ierror .gt. 0) stop



!
!---------- compute the kriging coefficients for
!           the local interpolation
!
!-------- Interpolation level # 1  
!
      ijcount=5
      iptl1=nvectpos
      iptcv1=kkcv    
      

       do kk=1,ilevref  
      
      if(imark.ne.1) then                                  
      
!------- compute the kriging coefficients for the level #1        
          jfin=2   
          do ii=1,4
             u1l1(ii,kk)=dble(0.0)
             u1l2(ii,kk)=dble(0.0)
          end do         
  

      do ii=1,3
         do jj=1,4
            cov1a(ii,jj)=cov1(ii,jj,kk) 
         end do
      end do
      

      call krig(1,ijcount,jfin,3,3,cov1a,vcxl1(kk))
      
      do ii=1,4
         u1l1(ii,kk)=u1(ii)
      end do   
      write(88)(u1(ii),ii=1,4)
      
      if(itype.ne.5) then
         vcxl1(kk)=dsqrt(one-vcxl1(kk))
         write(88)vcxl1(kk)
         vcxl1(kk)=sd1*vcxl1(kk)   
      else
      vcxl1(kk)=dsqrt(vcxl1(kk))
      write(88)vcxl1(kk)   
      end if
      
         

!      write(*,*)'level # 1:'
!      write(*,*)(u1l1(k),k=1,4)
!      write(*,*)dsqrt(vcxl1)

!------- compute the kriging coefficients for the level #2 

      do ii=1,10
         cov2a(ii)=cov2(ii,kk)
      end do

      call coefl2(itype,cov2a,u1,vcxl2(kk))       

      do ii=1,4
         u1l2(ii,kk)=u1(ii)
      end do   
      write(88)(u1(ii),ii=1,4)
      
      if(itype.ne.5) then
         vcxl2(kk)=dsqrt(one-vcxl2(kk))
         write(88)vcxl2(kk)
         vcxl2(kk)=sd1*vcxl2(kk)   
      else
      vcxl2(kk)=dsqrt(vcxl2(kk))
      write(88)vcxl2(kk)         
      end if
      
      else
!      
!---- read the kriging coefficients for the refinement levels
!
!
!------ level # 1
!      
         read(88)(u1l1(ii,kk),ii=1,4)
         read(88)vcxl1(kk)
         if(itype.ne.5) vcxl1(kk)=sd1*vcxl1(kk)
!
!------- level # 2
!
         read(88)(u1l2(ii,kk),ii=1,4)
         read(88)vcxl2(kk)
         if(itype.ne.5) vcxl2(kk)=sd1*vcxl2(kk)
      end if
      
                        
      end do
      close(88)

      write(*,*)

 



endsubroutine



