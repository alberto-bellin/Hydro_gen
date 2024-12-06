 subroutine gen_field2d
 use vars_common_hydrogen
 use mod_comb
 USE MOD_SEED_FROM_URANDOM
 
 
 implicit none

 real*8,allocatable,dimension(:)              :: r
 real*8,allocatable,dimension(:,:)            :: cond 
 integer (kind=4), allocatable                :: seed(:)

 integer  :: i,j,k,ndim_seed
 integer  :: igrid1
 integer  :: iiniz
 integer  :: ij,ij3
 integer  :: ip1,ip1s1,ip1s2,ipnod1,ipnod1s1,ipnod1s2
 integer  :: ips3
 integer  :: irv
 integer  :: ivectpos
 integer  :: isr
 integer  :: jfin,jiniz
 integer  :: n1l1,n2l1,n1finl2,n2finl2,jump1,jump2,inl1,inl2,infl2,kk
 integer  :: i1,j1,i0,j0,i2,j2,ijcount,i1a,j1a
 integer  :: dimu1
 integer  :: lag_x,lag_y



 real*8, parameter:: zero=0.0d0
 
 real*8   :: sumy,varlogy,condm,condv,xprint,yprint
 real*8   :: mean_r, var_r
 
 
 character*12 :: fileout
 
 complex ( kind = 8 ) c8_normal_01
 
 !real             :: gennor
 
 



 !----- pointers -------
 !
 !      zone       krig. coef.    cond. variance
 ! ----------------------------------------------
 ! |   side #1 |     ip1s1     |    ipnod1s1    |
 ! |   side #2 |     ip1s2     |    ipnod1s2    |
 ! |   side #3 |     ips3      |    iptcv3      |
 ! |   corners |     ip1       |    ipnod1      |
 ! |larger area|     iptkrs    |     ---        |
 !-----------------------------------------------
 !
 !
 !
 !----- generate the set of normally distributed random numbers---
 !

 allocate (cond(n1tot+1,n2tot+1))
 allocate (r(ngen))

 
 if(nseed.eq.0) then

    call random_seed(size = ndim_seed)
    allocate(seed(ndim_seed))
    call random_seed(get=seed)
    nseed = seed(1)

!    call SEED_FROM_URANDOM()
!    call random_seed(size = ndim_seed)
!    ndim_seed = 1
!    allocate(seed(ndim_seed))
!    call random_seed(get=seed)
!    nseed = seed(1)
    deallocate(seed)
    write(*,*)'generated seed:',nseed
 
 end if
 
 do k=1,ngen
       r(k)= c8_normal_01 ( nseed )
 end do


 
 !
 !--------------Initialize  RF to zero.--------
 !
 !
 do j=1,n2tot
    do i=1,n1tot
        cond(i,j)=zero
    end do
 end do  
 !
 !----- coarse grid ----
 !
 irv=1
 ivectpos=0
 kkcv=0  
 !
 !
 if(itype.eq.5) then !itype5if    
        !----- Fractal Field
        i=1
        do j=1,n2
           if(j.eq.1) then
              cond(1,1)=r(irv)            
           else
              kkcv=kkcv+1
              dimu1=nvectpos-ivectpos
              call comb(idiv,i,j,1,1,j,irv,dimu1,vcxsid(kkcv),u1sid(ivectpos+1:nvectpos),r,cond)    
           ivectpos=ivectpos+j-1
           end if
        end do
        !
        do i=2,n1
           do j=1,n2
             kkcv=kkcv+1
             dimu1=nvectpos-ivectpos
             call comb(idiv,i,j,1,1,n2,irv,dimu1,vcxsid(kkcv),u1sid(ivectpos+1:nvectpos),r,cond)              
             ivectpos=ivectpos+(i-1)*n2+j-1
           end do
        end do   
        !
  else  !itype5if-else    
        !      
        !----- regular random field
        !
        !
        !---- set the pointers for the kriging ceofficients and the conditional
        !     variances at the filed corners
        ip1=iptkrc
        ipnod1=iptcvc
        ips3=iptkr3
        !
        !------cycle through each row and column of the field-------
        !
        do 120 i=1,n1
           ip1s1=0
           ip1s2=iptkr2
           ipnod1s1=0
           ipnod1s2=iptcv2
           !
           do 121 j=1,n2
              !if we are at the first node, then we have no kriging to do,
              !create a random number.
              if(i.eq.1.and.j.eq.1) then
                 cond(1,1)=sd1*r(irv)
                 go to 121
              end if
              !
              IINIZ=I-NN1
              JINIZ=J-NN2
              jfin=j+nn2a
              !
              if(iiniz.lt.1) then
                 !
                 !that is, are we along side #3?
                  !
                  if(jiniz.ge.1.and.jfin.le.n2) then
                     !
                     !are we along side #3 but not on the corner?
                     !
                     !----- compute the velocities along the side # 3
                     !jiniz
                     !      |
                     !      v
                     !      --------------------* (i,j)       Note:  ijcount is the number of
                     !      |   |   |   |   |   |                    kriging coeffs for
                     !      |---|---|---|---|---|---|---| _          this position on
                     !      |   |   |   |   |   |   |   | ^          side #3
                     ! -----|---|---|---|---|---|---|---| nn1-----------
                     !      |   |   |   |   |   |   |   | v     edge of random field.
                     !   -> ----------------------------- -
                     ! iiniz                            ^
                     !      |<-----nn2--------->| nn2a  | jfin
                     !
                     iiniz=1
                     !loop through all kriging points to extract kriging coefficients
                     !
                     ij=(jfin-jiniz+1)*(i-1)+nn2
                     ij3=ij
                     !
                     dimu1=nvectpos-ips3            
                     call comb(idiv,i,j,iiniz,jiniz,jfin,irv,dimu1,vcxsid(i+iptcv3),u1sid(ips3+1:nvectpos),r,cond)
                     !  
                  else
                     !we are on a corner.
                     !
                     iiniz=1
                     jiniz=max(jiniz,1)
                     jfin=min(jfin,n2)
                     ij=(i-iiniz)*(jfin-jiniz+1)+j-jiniz
                     !
                     ipnod1=ipnod1+1
                     !            sc1=dsqrt(sigy-sigy*vcxkr(ipnod1))
                     !loop over all kriging coefficients for the corner to extract
                     
                     dimu1=nvectpos-ip1                     
                     call comb(idiv,i,j,iiniz,jiniz,jfin,irv,dimu1,vcxsid(ipnod1),u1sid(ip1+1:nvectpos),r,cond)
                     ip1=ip1+ij
                  endif
              else
                  !
                  !we are one side #1
                  !
                  if(jiniz.lt.1.and.jfin.le.n2) then
                    jiniz=1
                    ij=(jfin-jiniz+1)*nn1+j-1
                    !
                    ipnod1s1=ipnod1s1+1
!                   sc1=dsqrt(sigy-sigy*vcxsid1(ipnod1s1))
                    dimu1=nvectpos-ip1s1
                    call comb(idiv,i,j,iiniz,jiniz,jfin,irv,dimu1,vcxsid(ipnod1s1),u1sid(ip1s1+1:nvectpos),r,cond)
                    ip1s1=ip1s1+ij
                    go to 121
                  endif
                  !
                  !  we are on side #2
                  !
                  if(jfin.gt.n2) then
                    jfin=n2
                    ij=(jfin-jiniz+1)*nn1+j-jiniz
                    ipnod1s2=ipnod1s2+1
!                   sc1=dsqrt(sigy-sigy*vcxsid2(ipnod1s2))
                    dimu1=nvectpos-ip1s2
                    call comb(idiv,i,j,iiniz,jiniz,jfin,irv,dimu1,vcxsid(ipnod1s2),u1sid(ip1s2+1:nvectpos),r,cond)
                    ip1s2=ip1s2+ij
!                   nvectpos2=nvectpos2+(jfin-jiniz+1)*nn1+j-jiniz
!                   write(*,*)'side #2:',i,j,nvectpos2,ipnod1s2
                    go to  121
                   endif
                  !
                  ! we are in the greater kriging area
                  !
                  if(jiniz.ge.1.and.jfin.le.n2) then
                    ij=ijmax-1
                    dimu1=nvectpos-iptkrs
                    call comb(idiv,i,j,iiniz,jiniz,jfin,irv,dimu1,vcxsid(iptcvst),u1sid(iptkrs+1:nvectpos),r,cond)
!                   write(*,*)'steady:',i,j,ijcount
                  endif
              endif
           !
           ! Use extracted kriging coefficients to determine value of conductivity
           !
           121    continue
        ips3=ips3+ij3
        120   continue
        !
  endif!itype5-endif
  
       if(ilevref.gt.0) then  
 
   !
   !---- multistage grid refinement----
   !      
   !
   !---- modifications: August 22, 1997
   !
       n1l1=n1tot-idiv+1
       n2l1=n2tot-idiv+1  
       n1finl2=n1tot
       n2finl2=n2tot
       
   !
   !------------------------------
   !
       jump1=idiv  
       jump2=idiv/2
       inl1=1 
       inl2=idiv/2 +1
       infl2=inl2-jump2
!       
!       
       do kk=1,ilevref  
  
 
   !
   !------- first level refinement -----
 
 
        do 150 i1=inl1,n1l1,jump1
        i0=i1+jump1/2
  
        do 150 j1=inl1,n2l1,jump1   
        j0=j1+jump1/2  
 
        
        ijcount=0
        
        do i2=i1,i1+jump1,jump1
 
        do j2=j1,j1+jump1,jump1
        ijcount=ijcount+1 
 
        cond(i0,j0)=cond(i0,j0)+u1l1(ijcount,kk)*cond(i2,j2)
        end do
        end do 
 
   !
   !  Add the fluctuation
   !
        irv=irv+1
        cond(i0,j0)=cond(i0,j0)+vcxl1(kk)*r(irv)
 150    continue
        
 
   !       
   !---- second level refinement ----
   !      
        
        do 151 i1=inl2,n1finl2,jump1
        do 151 j1=inl2,n2finl2,jump1  
        
 
        i1a=i1
        j1a=j1-jump2
        if(j1a.ne.infl2) then   
  
 
               
        cond(i1a,j1a) =  cond(i1a,j1a)+cond(i1a,j1a-jump2)*u1l2(1,kk)+&
                        cond(i1a,j1a+jump2)*u1l2(3,kk)+&
                        cond(i1a-jump2,j1a)*u1l2(2,kk)+&
                        cond(i1a+jump2,j1a)*u1l2(4,kk)
   !
   ! Add the fluctuation
   !
        irv=irv+1
        cond(i1a,j1a)=cond(i1a,j1a)+vcxl2(kk)*r(irv)
        end if
 
        j1a=j1
        i1a=i1-jump2
 
        if(i1a.ne.infl2) then  
! 
! 
        cond(i1a,j1a) =  cond(i1a,j1a)+&
                        cond(i1a,j1a-jump2)*u1l2(1,kk)+&
                        cond(i1a,j1a+jump2)*u1l2(3,kk)+&
                        cond(i1a-jump2,j1a)*u1l2(2,kk)+&
                        cond(i1a+jump2,j1a)*u1l2(4,kk)
   !
   !  Add the fluctuation
   !
        irv=irv+1
        cond(i1a,j1a)=cond(i1a,j1a)+vcxl2(kk)*r(irv)
        end if
 151    continue 
  
         jump1=jump1/2
        inl1=inl1+jump1 
        jump2=jump2/2
        infl2=inl2
        inl2=inl2+jump2    
        n1finl2=n1finl2-jump2
        n2finl2=n2finl2-jump2 
        
        end do   
     
        if(irv.gt.ngen) then  
            write(*,*)'ngen=',ngen
            write(*,*)'code error: increase to',irv,'the parameter ngen'
            stop
        end if
 
       end if
!   !
!   !----- add the mean ----
!   !
! 
       sumy=dble(0.0)
       varlogy=dble(0.0)
       condm=dble(0.0)
       condv=dble(0.0)

  !      
  !*******************************************************
  !
  !        STATISTICS
  !
  !*******************************************************
  !
  !---------- spatial statistics------------
       
       do i=n1in,n1fin
          do j=n2in,n2fin
             condm=condm+cond(i,j)
             condv=condv+cond(i,j)*cond(i,j)
             cond(i,j)=cond(i,j)+cond10
          end do
       end do

       
       condm=condm/dble(nxy)
       condv=condv/dble(nxy-1) - condm*condm 
       condm=condm+cond10  
       
       if(itype.eq.5.and.answ.eq.'Y') then
       
       do i=n1in,n1fin
          do j=n2in,n2fin
             cond(i,j)=cond(i,j)-condm+fixcond
          end do
       end do                                 
       condm=fixcond
       end if
                
       
              write(20,9999)iNP,condm,condv
 9999         format(3x,i5,7x,2f15.8)      
  !
  !------ print the field -----
!  !
!  !
   fileout='real0000.dat'
 
 333   format(I1)
 334   format(I2)
 335   format(I3)
 336   format(I4)
! 	
! 		
 	if(iNP.lt.10) then
 	   write(fileout(8:8),333)iNP
 	end if
 
 	if(iNP.ge.10.and.iNP.lt.100) then
 	   write(fileout(7:8),334)iNP
 	end if
! 
 	if(iNP.ge.100.and.iNP.lt.1000) then
 	   write(fileout(6:8),335)iNP
 	end if  
! 
 	if(iNP.ge.1000.and.iNP.lt.10000) then
 	   write(fileout(5:8),336)iNP
 	end if  
! 
!! SEB commentato
!!       open(18,file=fileout)
 !! SEB
 
 !! SEB decommentato
   write(18,*)'replicate #',iNP
 ! SEB
             if(iformat.eq.1) then   
                   do i=n1in,n1fin
                    do j=n2in,n2fin
                        xprint=(j-1)*dx
                        yprint=(i-1)*dy
                        write(18,*)xprint,yprint,cond(i,j)
                    end do
                  end do
             else
                  do i=n1in,n1fin
                        write(18,866)(cond(i,j),j=n2in,n2fin)
                  end do
              end if
 
 
 866        format(6f20.10)
           close(18,status='keep')
!! SEB aggiunto
         close(20,status='keep')
!! SEB

!  compute the statistics of the field

 if(flag_stat.eq.1) then
     write(*,*)'compute the statistics of the field.....'

    lag_x=(n2fin-n2in)/2
    lag_y=(n1fin-n1in)/2
    call cov_spaz(n1fin-n1in+1,n2fin-n2in+1,lag_x,lag_y,cond(n1in:n1fin,n2in:n2fin))
 end if
 
 end subroutine
