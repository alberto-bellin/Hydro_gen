!c  file comb.f
!c
!c
!c-------------------------------------------------------------------
!c Author: Alberto Bellin  
!c Universita' di Trento,
! Dipartimento di Ingegneria Civile ed Ambientale,
! 38050-I Mesiano di Povo, TRENTO

!c e-mail: alberto.bellin@unitn.it 


!c Summary: It performs the linear combinations of selected field values
!c          to compute the conditional mean and generates the random
!c          deviate. 


!c Package Version: 2.0  April 1997 

!c Copyright: Alberto Bellin and Yoram Rubin.
!c
!c please cite the following paper in papers or reports that use
!c the present software:
!c Bellin A., Y. Rubin, Hydro_gen: A new random field generator
!c for correlated properties, Stochastic Hydrology and Hydraulics,
!c 10(4), 1996.
!c
!c
!c Permission is  granted to anyone to use and modify this packages provided
!c that:
!c i) the author is acknowledged;
!c ii) the use in any kind of research or job will be cited in the relative
!c papers or reports;
!c iii) anyone who uses the packages is under his/her own responsibility.
!c NO WARRANTY is given it is free of bugs and errors.
!c iv) The use or distribution must be free of charge.

!c Bug reports and hints are welcomed to the above e-mail address.
!c----------------------------------------------------------------------------

      module mod_comb
      contains
 
      subroutine comb(idiv,i,j,iiniz,jiniz,jfin,irv,dimu1,sc1,u1,r,cond)



 use vars_common_hydrogen,only:ngen,n1tot,n2tot

!      implicit double precision(a-g,o-z)

       implicit none


      integer,intent(in)  :: idiv
      integer,intent(in)  :: i
      integer,intent(in)  :: j
      integer,intent(in)  :: iiniz
      integer,intent(in)  :: jiniz
      integer,intent(in)  :: jfin
      integer,intent(in)  :: dimu1      


      real*8,intent(in)   :: sc1
      real*8,intent(in),dimension(ngen)      :: r
      real*8,intent(in),dimension(dimu1)      :: u1

      real*8,intent(out),dimension(n1tot+1,n2tot+1)   :: cond
      

      

       
       integer :: irv   !counter for random number extraction from the vector r
       integer :: ii0,jj0,iic,ii,j1,j1c,k,jjfin

!      dimension cond(idm,*),u1(*),r(*)  
      
             ii0=idiv*i-idiv+1
             jj0=idiv*j-idiv+1 
 
      
             k=0
              jjfin=jfin
             do  ii=iiniz,i
               if(ii.eq.i) jjfin=j-1
               do  j1=jiniz,jjfin
                 k=k+1
                 iic=idiv*ii-idiv+1
                 j1c=idiv*j1-idiv+1 

 
                 cond(ii0,jj0)=cond(ii0,jj0)+cond(iic,j1c)*u1(k)
               end do
             end do

              irv=irv+1
             cond(ii0,jj0)=cond(ii0,jj0)+sc1*r(irv)
      return
      end subroutine
      end module
