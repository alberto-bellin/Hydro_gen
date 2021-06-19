module vars_common_hydrogen
!---------------------------------------------------------------------------------------------
!  input
!---------------------------------------------------------------------------------------------
integer (kind =4) :: nseed      ! seed to initialize random generation
!integer (kind =4) :: nseed1    ! second seed
integer :: np                   ! number of Monte carlo realizations
integer :: imark                ! imark=1 to read the kriging coefficients from a file
integer :: itype                ! indicate the type of covariance function (0-5)
integer :: ilevref              ! number of refinement levels
integer :: iformat              ! define the  format of the ouput file containing the generated fields
integer :: ngen                 ! total # of nodes
integer :: n1tot                ! number of grid nodes along the y direction
integer :: n2tot                ! number of grid nodes along the x direction
integer :: ijmax                ! maximum dimension of the kriging system
integer :: flag_stat            ! flag to compute the statistics of the field (flag_stat=1)
integer, dimension (2)  :: n    ! define the number of grid nodes in both directions
!---------------------------------------------------------------------------------------------
!  input
!---------------------------------------------------------------------------------------------
real*8  :: dx,dy     ! final grid spacings 
real*8  :: ddx,ddy   ! grid spacing of the coarse grid
real*8  :: lx,ly     ! field dimensions
real*8  :: sigy      ! variance of the field
real*8  :: fixcond   ! mean value of the SRF imposed in each realization (read if answ=Y)
real*8  :: cond10    ! theorical mean value of the SRF (read if answ is not Y)
real*8  :: sclx,scly ! integral scales along the two principal directions
real*8  :: xsp       ! horizontal extension of the search neighborhood (left of the generation point)
real*8  :: ysp       ! vertical extension of the search neighnorhood
real*8  :: xspa      ! horizontal extension of the search neighborhood (rigth of the generation point
real*8  :: beta      ! exponent of the power law semivariogram
real*8  :: c         ! costant of the power law semivariogram
!---------------------------------------------------------------------------------------------
!  input
!---------------------------------------------------------------------------------------------
character*1 :: answ  !answ=Y to rescale the field 
character*30 :: file1 ! name output file (for the generated fields)
character*30 :: filecov
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
integer                           :: nvectpos
integer                           :: kkcv
integer                           :: iptcvst
integer                           :: n1in,n1fin,n2in,n2fin
integer                           :: iptkrc
integer                           :: iptcvc
integer                           :: iptkrs
integer                           :: iptkr2
integer                           :: iptcv2
integer                           :: iptkr3
integer                           :: iptcv3
integer                           :: iptl1
integer                           :: iptcv1

integer                           :: n1
integer                           :: n2
integer                           :: nn1
integer                           :: nn2
integer                           :: nn2a
integer                           :: nxy

integer                           :: idiv

integer                           :: iNP
!---------------------------------------------------------------------------------------------
real*8                            :: sd1
real*8                            :: vcx
!---------------------------------------------------------------------------------------------
real*8,allocatable,dimension(:)     :: u1
real*8,allocatable,dimension(:)     :: u1sid
real*8,allocatable,dimension(:)     :: vcxsid 
real*8,allocatable,dimension(:,:)   :: u1l1
real*8,allocatable,dimension(:,:)   :: u1l2
real*8,allocatable,dimension(:)     :: vcxl1
real*8,allocatable,dimension(:)     :: vcxl2
!real*8,allocatable,dimension(:)     :: r
!real*8,allocatable,dimension(:,:)   :: cond 
real*8,allocatable,dimension(:)     :: cov2a
real*8,allocatable,dimension(:,:,:) :: cov1 
real*8,allocatable,dimension(:,:)   :: cov2 
real*8,allocatable,dimension(:,:)   :: cov 
real*8,allocatable,dimension(:,:)   :: cov1a 

character(len=30),dimension(5)               :: covartype

end module
