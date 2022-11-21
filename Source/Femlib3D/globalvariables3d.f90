      module globalvariables3d
      use femtypes
!
!------------------------------------------------------------------------------
!    $Revision: 1.27 $
!    $Date: 2015/11/11 17:30:59 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  The globalvariables3d module defines frequently used variables to be avail-
!  able globally. Some variables only have a description without following
!  definition. If these will be used in a function/subroutine, the variable
!  names defined here should be used in the corresponding function/subroutine.
!
!------------------------------------------------------------------------------
!  GEOMETRY INFORMATION
!
!  Key-points:
!  -----------
!
!  Branches:
!  ---------
!
!  Regions:
!  --------
!     numdom    number of domains / regions
!     dommat    index of the material of a domain / region  allocate as: dommat(numdom)
!     domnames  Names of domains 
      integer (I4B) :: numdom
      integer (I4B), pointer :: dommat(:)=>null()
      character (len=40), pointer :: domnames(:)=>null()
      
!
!  Extents of geometry:
!  --------------------
!  This information cannot be accessed from NETGEN's geometry but must be ob-
!  tained by evaluation of the nodes in array nod.
!     xmin      minval(nod(1,:))
!     xmax      maxval(nod(1,:))
!     ymin      minval(nod(2,:))
!     ymax      maxval(nod(2,:))
!     zmin      minval(nod(3,:))
!     zmax      maxval(nod(3,:))
      real (DP) xmin, xmax, ymin, ymax, zmin, zmax
!
!------------------------------------------------------------------------------
!  MESH INFORMATION
!  (remark: not all variables used here are defined in globalvariables3D module
!           but should be used in common)
!
!  Numbers:
!  --------
!     nume      number of edges
!     numf      number of faces
!     numn      number of nodes
!     nums      number of surface elements
!     numv      number of volume elements
      integer (I4B), save :: nume, numn, nums, numv, numf
!
!  Nodes:
!  ------
!     nod       coordinates of nodes       allocate(nod(3,numn))
!               nod(1,:) => x-coordinates of nodes
!               nod(2,:) => y-coordinates of nodes
!               nod(3,:) => z-coordinates of nodes
      real (DP), pointer :: nod(:,:)=>null()
!
!  Mappings for volume elements (arrays):
!  --------------------------------------
!     ve        volume element -> edges    allocate(ve(6,numv))
!     vf        volume element -> faces    allocate(vf(4,numv))
!     vn        volume element -> nodes    allocate(vn(4,numv))
!     vv        neighbour information      allocate(vv(4,numv))
!               vv(i,j) is the element adjacent to j at element face i
!               (opposite to the vertex i). If element j has no neighbour at
!               face i,  vv(i,j)=0.
      integer (I4B), pointer :: ve(:,:)=>null(), vf(:,:)=>null(), vn(:,:)=>null(), vv(:,:)=>null()
!
!  Mappings for faces/edges to nodes (arrays):
!  -------------------------------------------
!     en        edge -> nodes              allocate(en(2,nume))
!     fn        face -> nodes              allocate(fn(3,numf))
!     sn        surface element -> nodes   allocate(sn(3,nums))
      integer (I4B), pointer :: en(:,:)=>null(), sn(:,:)=>null()
!
!  All mappings except for vv may be reverted by using the subroutine reversemap
!  leading to an array of pointers of type(ARRPTRI).
!
!  Number of natures
      integer (I4B) :: nnat
      
!  Boundary Conditions
! .--------------------------
!      numbc                number of boundary conditions
!      bcnames(:)           Text-specifier of boundary conditions                  bcnames(numbc)
!      bctype(:,:)          BC-type e.g. 0, 100, 200 ,...                          bctype(numbc,nnat)
!                             0 -  99: Dirichlet
!                           200 - 299: Neumannn
!                           300 - 399: Contour
!                           400 - 499: invisible
!      pvalue(:,:)          p-value of boundary condition                          pvalue(numbc,nnat)
!      qvalue(:,:)          q-value of boundary condition                          qvalue(numbc,nnat)
!      pvalue_vec(:,:,:)    p-value of vector boundary condition                   pvalue_vec(3,numbc,nnat)
!      qvalue_vec(:,:,:)    q-value of vector boundary condition                   qvalue_vec(3,numbc,nnat)
!      sbc(:)               BC-index of a surface element                          sbc(nums)
!      edgebc(:,)           BC-index of edges on surface                           edgebc(nume,nnat)
!      sfvertbc(:,:)        BC-index of vertices on surfaces                       sfvertbc(numn,nnat)

      integer (I4B) :: numbc
      integer (I4B), pointer :: sbc(:)=>null()
      integer (I4B), pointer :: bctype(:,:)=>null(), edgebc(:,:)=>null(), sfvertbc(:,:)=>null()
      complex (DPC), pointer :: pvalue(:,:)=>null(), qvalue(:,:)=>null()
      complex (DPC), pointer :: pvalue_vec(:,:,:)=>null(), qvalue_vec(:,:,:)=>null()
      character (len=40), pointer :: bcnames(:)=>null()
!
!  Auxillary mesh information:
!  ---------------------------
!     dom       assignment of element to domain / region
!     eltype    type of element/ shape function used
!     bcs is for each nature bcs(:,nnat)
      integer (I4B), pointer :: bcs(:,:)=>null()
      integer (I4B), pointer :: dom(:)=>null()
      character (len=16) :: eltype = 'NEDELEC'

!
!------------------------------------------------------------------------------
!  SOLUTION INFORMATION
!  (remark: not all variables used here are defined in globalvariables3D module
!           but should be used in common)
!
!  Degrees of freedom:
!  -------------------
!     numdof    total number of DOF
!     edof      DOFs of edge
!     egdof     edge: local DOF -> global DOF
!     fdof      DOFs of face
!     fgdof     face: local DOF -> global DOF
!     vdof      DOFs of volume element
!     vgdof     volume element: local DOF -> global DOF     allocate(vgdof(numv,nnat))
      integer (I4B), save :: numdof
      integer (I4B) , pointer :: natdof(:)=>null()
      type (ARRPTRI), pointer :: vgdof(:,:)=>null()
!
!  Polynomial degree:
!  -------------------
!
!     fp        polynomial degree of face              allocate(fp(numf,nnat))
!               for each nature
!     vp        polynomial degree of volume element    allocate(vp(numv,nnat))
!               for each nature
!     ep        polynomial degree of edges             allocate(ep(nume,nnat))
!               for each nature (if edge is on a boundary, the ep(i,nat=1) is negative)


      integer (I4B), pointer :: vp_temp(:,:)=>null(), vp(:,:)=>null(), fp(:,:)=>null(), ep(:,:)=>null()
!
!  Solution vector:
!  ----------------
!     x         solution vector with length numdof     allocate(x(numdof))
!     resn      average residual at nodes              allocate(resdof(numdof))
      complex (DPC), pointer :: x(:)=>null()
      real (DP), pointer :: resdof(:)=>null()
      
!  Variables for an auxiliary solution (corresponding to definitions above):
!  ----------------
!     x_aux     auxiliary solution vector with length numdof_aux     allocate(x_aux(numdof_aux))
!     nnat_aux
!     numdof_aux
!     vgdof_aux
!     vp_aux
!
      complex (DPC), pointer    :: x_aux(:)=>null()
      integer (I4B), pointer    :: vp_aux(:,:)=>null()
      integer (I4B), save       :: numdof_aux
      type (ARRPTRI), pointer   :: vgdof_aux(:,:)=>null()
      integer (I4B)             :: nnat_aux
!
!------------------------------------------------------------------------------
!  OTHER GLOBAL SETTINGS AND PARAMETERS
!
!  Auxillary:
!  ----------
!     omega     angular frequency
!     physics   determines PDE type, only wave optics atm. (TEWAVE,TMWAVE)
!     whatsystem   default value for operating system
      integer (I4B) :: arg_num_solver3d
      real (DP) :: omega
      character (len=10) , allocatable :: arg_solver3d(:)
      character (len=7)  :: whatsystem = 'Windows' 
      character (len=20) :: physics

!  Parameters:
!  -----------
!     polymaxsc max. available order of scalar (Lagrange)shape functions
!     polymaxc  max. available order of vectorial (Nédélec) shape functions
!     c0        free space speed of light
!     cs0       speed sound in air at 20 deg. C
!     eps0      free space permittivity
!     mu0       free space permeability
      integer (I4B), parameter :: polymaxsc = 9, polymaxv=5
      real (DP), parameter :: cs0 = 340.5_DP
      real (DP), parameter :: pi = 3.141592653589793238462643_DP
      real (DP), parameter :: Kboltz = 1.3806488E-23_DP                ! Boltzmann's constant
      real (DP), parameter :: KboltzEv = 8.617343e-5_DP                ! Boltzmann's constant (in eV/K)
      real (DP), parameter :: mu0 = 0.1256637061435917295385057e-5_DP  ! free space permeability
      real (DP), parameter :: c0 = 299792458._DP                       ! speed of light
      real (DP), parameter :: eps0 = .8854187817620389850536565e-11_DP ! free space permittivity
      real (DP), parameter :: gravc = 6.67428e-11_DP                   ! gravitational constant
      real (DP), parameter :: Elch = 1.60217646e-19_DP                 ! elementary charge (in Coulombs)
      real (DP), parameter :: Avogadro = 6.02214e23_DP                 ! The number of molecules in a mole (in 1/mol)
!
!------------------------------------------------------------------------------
      end module globalvariables3d
