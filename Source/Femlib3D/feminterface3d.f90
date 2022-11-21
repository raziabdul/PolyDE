module feminterface3d
!------------------------------------------------------------------------------
!
!    $Revision: 1.77 $
!    $Date: 2015/11/11 17:36:23 $
!    $Author: m_kasper $
!
!------------------------------------------------------------------------------
interface
      subroutine adaptation3D (epsgl,errcode)
      use femtypes
      implicit none
      integer (I4B) :: errcode
      real (DP)     :: epsgl
      intent (out)  :: epsgl, errcode
      end subroutine adaptation3D
end interface


interface
      subroutine apply_dirichlet_bc(elem, vert, nff, nffsum, polylo, polyhi, a, b)
      use femtypes
      implicit none
      integer (I4B) :: elem, nff(:), nffsum(:), polylo(:), polyhi(:)
      real (DP) :: vert(3,4)
      complex (DPC) :: a(:,:), b(:)
      intent (in) :: elem, vert, nff, nffsum, polylo, polyhi
      intent (inout) :: a, b
      end subroutine apply_dirichlet_bc
end interface


interface
      subroutine apply_neumann_bc(elem, vert, nff, nffsum, polylo, polyhi, a, b)
      use femtypes
      implicit none
      integer (I4B) :: elem, nff(:), nffsum(:), polylo(:), polyhi(:)
      real (DP) :: vert(3,4)
      complex (DPC) :: a(:,:), b(:)
      intent (in) :: elem, vert, nff, nffsum, polylo, polyhi
      intent (inout) :: a, b
      end subroutine apply_neumann_bc
end interface


interface
      pure function area3d(vert)
      use femtypes
      implicit none
      real (DP) :: area3d, vert(3,3)
      intent(in) :: vert
      end function area3d
end interface


interface
      subroutine assembly3D(lower,upper,diag,acsr,rhs,ia,ja,jacobi,csr, &
     &                      matvar,symmetric)
      use femtypes
      implicit none
      integer (I4B), pointer :: ia(:), ja(:)
      complex (DPC), pointer :: lower(:), upper(:), diag(:), rhs(:), acsr(:)
      logical jacobi, csr, matvar, symmetric
      intent (in) ::  jacobi, csr, matvar, symmetric
      end subroutine assembly3D
end interface


interface
      subroutine bcstovv
      use femtypes
      implicit none
      end subroutine bcstovv
end interface


interface
      subroutine BisectHull(e_0,es,ev,e_len,v_mark,n_0_in)
      use femtypes
      implicit none
      integer  (I4B)               :: e_0
      integer  (I4B),  allocatable :: v_mark(:)
      real      (DP),      pointer :: e_len(:)
      type (ARRPTRI),      pointer :: ev(:), es(:)
      integer  (I4B),     optional :: n_0_in
      intent    (in)               :: e_0, n_0_in
      intent (inout)               :: ev, es, e_len, v_mark
      end subroutine BisectHull
end interface


interface
      subroutine BisectTET(t_0,es,ev,e_len,v_mark)
      use femtypes
      implicit none
      integer (I4B)               :: t_0
      integer (I4B),  allocatable :: v_mark(:)
      real     (DP),      pointer :: e_len(:)
      type(ARRPTRI),      pointer :: ev(:), es(:)
      intent   (in)               :: t_0
      intent(inout)               :: ev, es, e_len, v_mark
      end subroutine BisectTET
end interface


interface
      subroutine bcs2vert(ok)
      use femtypes
      implicit none
      logical ok
      intent (out) :: ok
      end subroutine
end interface


interface
      subroutine calctensor3D(tensor,vect1,vect2)
      use femtypes
      implicit none
      real (DP) :: vect1(3), vect2(3)
      complex (DPC) :: tensor(3,3)
      intent(in) :: vect1, vect2
      intent(inout) :: tensor
      end subroutine calctensor3D
end interface


interface
      subroutine coloringsc3D(indx,offset)
      use femtypes
      implicit none
      integer (I4B), allocatable :: indx(:), offset(:)
      intent (out) :: indx, offset
      end subroutine coloringsc3D
end interface


interface
      elemental subroutine compareandswap(a,b)
      use femtypes
      implicit none
      integer (I4B) :: a, b
      intent (inout) :: a, b
      end subroutine compareandswap
end interface


interface
     subroutine CSR_residual(a,b,x,ia,ja,n,res)
        use feminterface, only:
        use femtypes
        implicit none
        complex (DPC) a(:), x(:), b(:), res(:)
        integer (I4B) n, ia(:), ja(:)
        intent (in) :: a, b, x, ia, ja, n
        intent (out) :: res
      end subroutine CSR_residual
end interface


interface
      subroutine curlshapev3D( l, gl, polylo, polyhi, vec, axs )
      Use femtypes
      implicit none
      integer (I4B) :: polylo, polyhi, axs
      real (DP)     :: l(4), gl(3,4)
      real (DP)     :: vec(:)
      intent (in)   :: l, gl, polylo, polyhi, axs
      intent (out)  :: vec
      end subroutine curlshapev3D
end interface


interface
      subroutine dgshapesc3D(l, gl, nu, polylo, polyhi, vec, errcode)
      use femtypes
      implicit none
      integer (I4B) polylo, polyhi, errcode
      real (DP) l(4), gl(3,4)
      complex (DPC) nu(:,:), vec(:)
      intent (in) :: l, gl, polylo, polyhi
      intent (out) :: vec, errcode
      end subroutine dgshapesc3D
end interface


interface
      pure function diangles(vert)
      use femtypes
      implicit none
      real   (DP) :: diangles(6), vert(3,4)
      intent (in) :: vert
      end function diangles
end interface


interface
      subroutine dirichlet3D( num, xyzs, dbc )
      use femtypes
      implicit none
      integer (I4B) :: num
      real (DP) :: xyzs(3)
      complex(DPC) :: dbc(3)
      intent (in)   :: num, xyzs
      intent (out)  :: dbc
      end subroutine dirichlet3D
end interface


interface
      function EdgeLength( n_0 , n_t )
      use femtypes                   
      implicit none                  
      integer (I4B)               :: n_0, n_t
      real     (DP)               :: EdgeLength
      intent   (in)               :: n_0, n_t
      end function EdgeLength
end interface


interface
      subroutine element_draw_vtk(ellist)
      use femtypes
      implicit none
      integer(I4B) :: ellist(:)
      intent (in) :: ellist
      end subroutine element_draw_vtk
end interface


interface
      subroutine elementmatrix3D(elem, vpelem, jacobi, full, matvar, a, b, nff, errcode)
      use femtypes
      implicit none
      integer (I4B) elem
      integer (I4B) vpelem, nff, errcode
      complex (DPC), allocatable :: a(:,:), b(:)
      logical jacobi, full, matvar
      intent (in) :: elem, vpelem, jacobi, full, matvar
      intent (out) :: nff, a, b
      end subroutine elementmatrix3D
end interface


interface
      subroutine elementmatrix_scalar3D(elem, jacobi, full, matvar, a, b, nff, nffsum, errcode)
      use femtypes
      implicit none
      integer (I4B) elem
      integer (I4B), allocatable :: nff(:), nffsum(:)
      integer (I4B) errcode
      complex (DPC), allocatable :: a(:,:), b(:)
      logical jacobi, full, matvar
      intent (in) :: elem, jacobi, full, matvar
      intent (out) :: a, b, nff, nffsum, errcode
      end subroutine elementmatrix_scalar3D
end interface


interface
      subroutine element_mark3D(crt,mark_confirm,element_mark,mark_depth,error_bound,marking_type,num_marked,produce_plot)
      use femtypes
      implicit none
      integer (I4B)               :: mark_depth, num_marked
      integer (I4B)               :: element_mark(:)
      integer (I4B)               :: mark_confirm(:)
      real     (DP)               :: error_bound
      real     (DP)               :: crt(:)
      character (len=*), optional :: produce_plot
      character (len=16)          :: marking_type
      intent   (in)               :: mark_depth, error_bound
      intent(inout)               :: element_mark, crt
      end subroutine element_mark3D
end interface


interface
      subroutine estimate_error3D(v_res,v_res_norm)
      use femtypes
      implicit none
      real (DP) :: v_res(:,:), v_res_norm(:)
      intent(out) :: v_res, v_res_norm
      end subroutine estimate_error3D
end interface


interface
    subroutine export2matlab(datatype,val)
    use femtypes
    implicit none
    real      (DP), optional    :: val(:,:)
    character (len=*)           :: datatype
    end subroutine export2matlab
end interface


interface
      subroutine exportdata3D(fieldtype,phi,origin,p1,p2,grid1,grid2)
      use femtypes
      implicit none
      integer (I4B) :: grid1, grid2
      real (DP) :: phi, origin(3), p1(3), p2(3)
      character (len=10) :: fieldtype
      intent (in) :: fieldtype, phi, origin, p1, p2, grid1, grid2
      end subroutine exportdata3D
end interface


interface
      subroutine exportvtk(fieldtype,phi)
      use femtypes
      implicit none
      real (DP) :: phi
      character (len=10) :: fieldtype
      intent (in) :: fieldtype, phi
      end subroutine exportvtk
end interface


interface
      subroutine exportvtkres(fieldtype,phi)
        use femtypes
        implicit none
        real (DP) :: phi
        character (len=10) :: fieldtype
        intent (in) :: fieldtype, phi
      end subroutine exportvtkres
end interface


interface
      subroutine exportvtkvar(fieldtype,phi)
      use femtypes
      implicit none
      real (DP) :: phi
      character (len=10) :: fieldtype
      intent (in) :: fieldtype, phi
      end subroutine exportvtkvar
end interface


interface
      subroutine field3D(elem,lambda,u,curlu)
      use femtypes
      implicit none
      integer (I4B) elem
      real (DP) lambda(4)
      complex (DPC) u(3), curlu(3)
      intent (in) ::  elem, lambda
      intent (out) ::  u, curlu
      end subroutine field3D
end interface


interface
      subroutine fieldsc3D(elem,lambda,typ,z,soln,u,alphau,gammau,gradu,betagradu,gammagradu,nugradu,dgradu,f,g)
      use femtypes
      implicit none
      integer (I4B) elem
      real (DP) lambda(4)
      logical typ(5)
      complex (DPC) z(:,:)
      complex (DPC), optional :: u(:), alphau(:,:), gammau(:,:,:), dgradu(:,:), f(:), g(:,:)
      complex (DPC), optional :: betagradu(:,:), gammagradu(:,:), nugradu(:,:,:)
      complex (DPC), optional :: gradu(:,:)
      complex (DPC), optional :: soln(:)
      intent (in) ::  elem, lambda, typ, soln
      intent (out) ::  z, u, alphau, gammau, gradu, betagradu, gammagradu, nugradu, dgradu, f, g
      end subroutine fieldsc3D
end interface


interface fieldsc3D_simple
      subroutine fieldsc3D_simple(elem,lambda,u,gu)
      use femtypes
      implicit none
      integer (I4B)           :: elem
      real (DP)               :: lambda(4)
      complex (DPC)           :: u(:)
      complex (DPC), optional :: gu(:,:)
      intent (in)             :: elem, lambda
      intent (out)            :: u
      end subroutine fieldsc3D_simple
end interface

interface fieldsc3D_simple_aux
      subroutine fieldsc3D_simple_aux(elem,lambda,u,gu)
      use femtypes
      implicit none
      integer (I4B)           :: elem
      real (DP)               :: lambda(4)
      complex (DPC)           :: u(:)
      complex (DPC), optional :: gu(:,:)
      intent (in)             :: elem, lambda
      intent (out)            :: u
      end subroutine fieldsc3D_simple_aux
end interface


interface fieldquantity3d
      subroutine fieldquantity3d_scalar(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
      use femtypes
      implicit none
      integer (I4B) :: elem
      real (DP) :: lambda(4), phi
      complex (DPC) :: zs
      character (len=*) :: fieldtype, unit, descriptor
      logical :: ok
      intent (in) :: elem, fieldtype, lambda, phi
      intent (out) :: descriptor, unit, zs, ok
      end subroutine fieldquantity3d_scalar

      subroutine fieldquantity3D_vector(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
      use femtypes
      implicit none
      integer (I4B) :: elem
      real (DP) :: lambda(4), phi
      complex (DPC) :: zs(:)
      character (len=*) :: fieldtype, unit, descriptor
      logical :: ok
      intent (in) :: elem, fieldtype, lambda, phi
      intent (out) :: descriptor, unit, zs, ok
      end subroutine fieldquantity3d_vector

      subroutine fieldquantity3d_tensor(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
      use femtypes
      implicit none
      integer (I4B) :: elem
      real (DP) :: lambda(4), phi
      complex (DPC) :: zs(:,:)
      character (len=*) :: fieldtype, unit, descriptor
      logical :: ok
      intent (in) :: elem, fieldtype, lambda, phi
      intent (out) :: descriptor, unit, zs, ok
      end subroutine fieldquantity3d_tensor
end interface


interface
      subroutine findelement3D(xyzt,nod,vn,vv,numv,ielem,found)
      use femtypes
      real (DP) :: xyzt(:), nod(:,:)
      integer (I4B) :: vn(:,:), vv(:,:), numv, ielem
      logical :: found
      intent (in) :: xyzt, nod, vn, vv, numv
      intent (out) :: ielem, found
      end subroutine findelement3D
end interface


interface
      subroutine findelement_x3D(xyzt,nod,vn,numv,ielem,found)
      use femtypes
      real (DP) :: xyzt(:), nod(:,:)
      integer (I4B) :: vn(:,:), numv, ielem
      logical :: found
      intent (in) :: xyzt, nod, vn, numv
      intent (out) :: ielem, found
      end subroutine findelement_x3D
end interface


interface
      subroutine get_modres(v_res,mod_res,vol_poly)
      use femtypes
      implicit none
      real     (DP)              :: v_res(:), mod_res(:)
      integer  (I4B)             :: vol_poly(:)
      intent(in) :: v_res, vol_poly
      intent(out) :: mod_res
      end subroutine get_modres
end interface


interface
      subroutine get_pconfirm(sing_val,max_vp,pconfirm)
      use femtypes
      implicit none
      real     (DP)                :: sing_val(:)
      integer (I4B)                :: max_vp, pconfirm(:)
      end subroutine get_pconfirm
end interface


interface
      pure subroutine get2Dintpolpoints(norder, lambda, npkt, errcode)
      use femtypes
      implicit none
      integer (I4B) norder,  npkt, errcode
      real (DP), allocatable :: lambda(:,:)
      intent (in) :: norder
      intent (out) :: errcode, npkt, lambda
      end subroutine get2Dintpolpoints
end interface


interface
      pure subroutine get3Dintegpoints(intorder, npkt, weig, lambda, errcode)
      use femtypes
      implicit none
      integer (I4B) intorder, npkt, errcode
      real (DP), allocatable :: weig(:), lambda(:,:)
      intent (in) :: intorder
      intent (out) :: npkt, errcode, weig, lambda
      end subroutine get3Dintegpoints
end interface


interface
      subroutine getbcval( face, nature, xyzs, pval, qval )
      use femtypes
      implicit none
      integer (I4B) :: face, nature
      real (DP) :: xyzs(3)
      complex(DPC) :: pval
      complex(DPC), optional :: qval(:)
      intent (in)   :: face, nature, xyzs
      intent (out)  :: pval, qval
      end subroutine getbcval
end interface


interface
      subroutine getbcval_vec( face, nature, xyzs, pval, qval )
      use femtypes
      implicit none
      integer (I4B) :: face, nature
      real (DP) :: xyzs(3)
      complex(DPC) :: pval(3)
      complex(DPC), optional :: qval(3,3)
      intent (in)   :: face, nature, xyzs
      intent (out)  :: pval, qval
      end subroutine getbcval_vec
end interface


interface
      subroutine getCriticalPath(t_0,CP,n_0,ev,e_len)
      use femtypes
      implicit none
      integer  (I4B)                    :: t_0
      integer  (I4B),           pointer :: CP(:), n_0(:)
      real     (DP),            pointer :: e_len(:)
      type(ARRPTRI),            pointer :: ev(:)
      intent    (in)                    :: t_0
      intent   (out)                    :: CP, n_0
      intent (inout)                    :: ev, e_len
      end subroutine getCriticalPath
end interface


interface
      subroutine getepbc
      use femtypes
      implicit none
      end subroutine getepbc
end interface


interface
      subroutine getes(es)
      use femtypes
      implicit none
      type(ARRPTRI),     pointer :: es(:)
      intent(out)                :: es
      end subroutine getes
end interface


interface
      subroutine getfe(fe)
      use femtypes
      implicit none
      integer  (I4B), allocatable :: fe(:,:)
      intent   (out)              :: fe
      end subroutine getfe
end interface


interface
      pure subroutine getgl(vert, gl, n, t1, t2)
      use femtypes
      implicit none
      real (DP) :: vert(3,4), gl(3,4)
      real (DP), optional :: n(3,4), t1(3,4), t2(3,4)
      intent (in) :: vert
      intent (out) :: gl, n, t1, t2
      end subroutine getgl
end interface


interface
      subroutine getkpv(mod_res,keyPointValues)
      use femtypes
      implicit none
      real     (DP)             :: mod_res(:), keyPointValues(:)
      end subroutine getkpv
end interface


interface
      subroutine getse(se)
      use femtypes
      implicit none
      integer  (I4B), allocatable :: se(:,:)
      intent   (out)              :: se
      end subroutine getse
end interface


interface
      subroutine get_starting_solution (x_0)
      use femtypes
      use matconstants
      implicit none
      complex (DPC), pointer    :: x_0(:)
      intent (inout)            :: x_0
      end subroutine get_starting_solution
end interface


interface
      subroutine getvf(numf,vf)
      use femtypes
      implicit none
      integer (I4B) :: numf
      integer (I4B), pointer :: vf(:,:)
      intent (out) :: numf, vf
      end subroutine getvf
end interface


interface
      subroutine gradshapesc3D( l, gl, polylo, polyhi, vec)
      use femtypes
      implicit none
      integer (I4B) :: polylo, polyhi
      real (DP)     :: l(4), gl(4)
      real (DP)     :: vec(:)
      intent (in)   :: l, gl, polylo, polyhi
      intent (out)  :: vec
      end subroutine gradshapesc3D
end interface


interface
      subroutine h_adapt3D(astep,v_res,ev,es,e_len)
      use femtypes
      implicit none
      integer (I4B)          :: astep
      real     (DP)          :: v_res(:,:)
      real     (DP), pointer :: e_len(:)
      type(ARRPTRI), pointer :: ev(:), es(:)
      intent(in) :: astep
      intent(out) :: v_res
      intent(inout) :: ev, es, e_len
      end subroutine h_adapt3D
end interface


interface
      subroutine hp_adapt3D(astep,v_res,ev,es,e_len)
      use femtypes
      implicit none
      integer (I4B)          :: astep
      real     (DP)          :: v_res(:,:)
      real     (DP), pointer :: e_len(:)
      type(ARRPTRI), pointer :: ev(:), es(:)
      intent(in) :: astep
      intent(out) :: v_res
      intent(inout) :: ev, es, e_len
      end subroutine hp_adapt3D
end interface


interface
      subroutine hp_alg_F8R(astep,v_res,ev,es,e_len)
      use femtypes
      implicit none
      integer (I4B)          :: astep
      real     (DP)          :: v_res(:,:)
      real     (DP), pointer :: e_len(:)
      type(ARRPTRI), pointer :: ev(:), es(:)
      intent(in) :: astep
      intent(inout) :: v_res, ev,es,e_len
      end subroutine hp_alg_F8R
end interface


interface
      subroutine hmesh_refine3D(v_mark,es,ev,e_len,num_refined)
      use femtypes
      implicit none
      integer (I4B)               :: num_refined
      integer (I4B),  allocatable :: v_mark(:)
      real     (DP),      pointer :: e_len(:)
      type(ARRPTRI),      pointer :: ev(:), es(:)
      intent(out)                 :: num_refined
      intent(inout)               :: es, ev, e_len, v_mark
      end subroutine hmesh_refine3D
end interface


interface
      subroutine initialize3D()
      use femtypes
      implicit none
      end subroutine initialize3D
end interface


interface
      subroutine lam2xyz(lambda, elem, xyzt, nod, vn)
      use femtypes
      implicit none
      integer (I4B) elem, vn(:,:)
      real (DP) xyzt(3), lambda(4), nod(:,:)
      intent (in) :: lambda, nod, elem, vn
      intent (out) :: xyzt
      end subroutine lam2xyz
end interface


interface
      subroutine lambcoef(v1, v2, v3, v4, a, b, c, d)
      use femtypes
      implicit none
      real (DP) v1(3), v2(3), v3(3), v4(3)
      real (DP) a(4), b(4), c(4), d(4)
      intent (in) :: v1, v2, v3, v4
      intent (out) :: a, b, c, d
      end subroutine lambcoef
end interface


interface
      subroutine LCSR_residual(lower,upper,diag,b,x,ia,ja,n,res)
      use femtypes
      implicit none
      complex (DPC) lower(:), upper(:), diag(:), x(:), b(:), res(:)
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: upper, lower, diag, b, x, ia, ja, n
      intent (out) :: res
      end subroutine LCSR_residual
end interface


interface
      subroutine meshquality3d(vn,numv,nod)
      use femtypes
      implicit none
      integer (I4B) :: vn(:,:), numv
      real (DP) :: nod(:,:)
      intent(in) :: vn, numv, nod
      end subroutine meshquality3d
end interface


interface
      subroutine meshqualityvtk()
      use femtypes
      implicit none
      end subroutine meshqualityvtk
end interface


interface
      subroutine nedelecdof(numf, vf)
      use femtypes
      implicit none
      integer (I4B), pointer :: vf(:,:)
      integer (I4B), intent (in) :: numf
      end subroutine nedelecdof
end interface


interface
      subroutine p_adapt3D(astep,v_res)
      use femtypes
      implicit none
      integer (I4B)          :: astep
      real     (DP)          :: v_res(:,:)
      intent (in) :: astep
      intent (out) :: v_res
      end subroutine p_adapt3D
end interface


interface
      subroutine petscNonlinear(xn)
      use femtypes
      implicit none
      complex (DPC), pointer ::  xn(:)
      intent (inout) :: xn
      end subroutine petscNonlinear
end interface

interface
      subroutine pmesh_refine3D(v_mark,nat,num_refined,ref_vp,pconfirm)
      use femtypes
      implicit none
      integer (I4B), optional              :: pconfirm(:,:)
      integer (I4B)                        :: nat, num_refined
      integer (I4B)                        :: v_mark(:)
      logical                              :: ref_vp
      end subroutine pmesh_refine3D
end interface


interface pdecoeff3D
! for vector-valued basis functions
      subroutine pdecoeff3D_vector(elem, xyzs, nu, alpha, beta, f1, f2)
      use femtypes
      implicit none
      integer (I4B) elem
      real(DP) :: xyzs(3)
      complex (DPC) :: nu(3,3), alpha(3,3), beta(3)
      complex (DPC) :: f1(3), f2(3)
      intent (in) :: elem, xyzs
      intent (out) :: nu, alpha, beta, f1, f2
      end subroutine pdecoeff3D_vector
!
! for scalar-valued basis functions
      subroutine pdecoeff3D_scalar(elem, xyzs, nu, alpha, beta, f1, f2, gamma)
      use femtypes
      implicit none
      integer (I4B) elem
      real(DP) :: xyzs(3)
      complex (DPC) :: nu(:,:,:,:), alpha(:,:), beta(:,:,:)
      complex (DPC) :: f1(:), f2(:,:)
      complex (DPC), optional :: gamma(:,:,:)
      intent (in) :: elem, xyzs
      intent (out) :: nu, alpha, beta, f1, f2, gamma
      end subroutine pdecoeff3D_scalar
end interface pdecoeff3D


interface
      subroutine plot_marking(plot_input)
      use femtypes
      implicit none
      real (DP)                   :: plot_input(:)
      intent(in)                  :: plot_input
      end subroutine plot_marking
end interface


interface
      subroutine pointvalue3D
      use femtypes
      implicit none
      end subroutine pointvalue3D
end interface


interface
      subroutine polyordervtk()
      use femtypes
      implicit none
      end subroutine polyordervtk
end interface


interface
      subroutine preassemb3D()
      use femtypes
      implicit none
      end subroutine preassemb3D
end interface


interface
      subroutine prepare_h_adapt3D(ev,es,e_len)
      use femtypes
      implicit none
      real     (DP),     pointer :: e_len(:)
      type(ARRPTRI),     pointer :: ev(:), es(:)
      intent(out)                :: ev, es, e_len
      end subroutine prepare_h_adapt3D
end interface


interface
      subroutine prepare_hp_adapt3D(ev,es,e_len)
      use femtypes
      implicit none
      real     (DP),     pointer :: e_len(:)
      type(ARRPTRI),     pointer :: ev(:), es(:)
      intent(out)                :: es, ev, e_len
      end subroutine prepare_hp_adapt3D
end interface


interface
      subroutine prepare_p_adapt3D
      use femtypes
      implicit none
      end subroutine prepare_p_adapt3D
end interface


interface
      subroutine principlevalues3D(vx,vy,vz,vxy,vxz,vyz,lambda1,lambda2,lambda3,angle)
      use femtypes
      implicit none
      real (DP) :: vx,vy,vz,vxy,vxz,vyz,lambda1,lambda2,lambda3
      real (DP), optional :: angle
      intent (in) :: vx,vy,vz,vxy,vxz,vyz
      intent (out) :: lambda1, lambda2, lambda3, angle
      end subroutine principlevalues3D
end interface


interface
      subroutine readbc(bcnames, numbc, nnat, hasnames, bctype, pvalue, qvalue, ok)
      use femtypes
      implicit none
      integer(I4B) :: numbc, nnat, bctype(:,:)
      complex(DPC) :: pvalue(:,:), qvalue(:,:)
      character (len=*) :: bcnames(:)
      logical :: ok, hasnames
      intent (in) :: numbc, nnat, hasnames
      intent (out) :: bctype, pvalue, qvalue, ok
      intent (inout) :: bcnames
      end subroutine readbc
end interface


interface
      subroutine readmat3D(matname,found,matindex)
      use femtypes
      implicit none
      character (len=*) matname
      integer (I4B) matindex
      logical found
      intent (out) :: found, matindex
      intent (in) :: matname
    end subroutine readmat3D
end interface


interface
      subroutine readng(meshfile, ok)
      use femtypes
      implicit none
      character(len=*) :: meshfile
      logical ok
      intent (in) :: meshfile
      intent (out) :: ok
      end subroutine readng
end interface


interface
      subroutine readunv(meshfile,ok)
      use femtypes
      implicit none
      character(len=*) :: meshfile
      logical ok
      intent (in) :: meshfile
      intent (out) :: ok
      end subroutine readunv
end interface


interface
      subroutine reversemap(column,array,arrptr,nmin,nmax)
      use femtypes
      implicit none
      integer (I4B) :: nmin, nmax
      integer (I4B) :: array(:,:)
      type (ARRPTRI), pointer :: arrptr(:)
      logical :: column
      intent(in) :: column
      intent(out) :: nmin, nmax
      end subroutine reversemap
end interface


interface
      subroutine residual3D(res,sumref)
      use femtypes
      implicit none
      real(DP), pointer :: res(:)
      real(DP), intent(out) :: sumref
      end subroutine residual3D
end interface


interface
      subroutine residualsc3D(ext,int,matvar,errcode,res,sumref,resext)
      use femtypes
      implicit none
      integer (I4B) errcode
      real(DP):: res(:,:), sumref(:)
      real(DP), pointer, optional :: resext(:,:,:)
      logical :: ext, int, matvar
      intent (in) :: ext, int, matvar 
      intent (out) :: errcode, res, sumref, resext
      end subroutine residualsc3D
end interface


interface
      subroutine setstandardvalues3D
      use femtypes
      implicit none
      end subroutine setstandardvalues3D
end interface


interface
      subroutine scalardof(numf, vf)
      use femtypes
      implicit none
      integer (I4B), pointer :: vf(:,:)
      integer (I4B), intent (in) :: numf
      end subroutine scalardof
end interface


interface
      subroutine shapefunctionv3D(lambda, vert, polylo, polyhi,   &
     &  nff, curl, xsi, cxsi, errcode)
      use femtypes
      implicit none
      integer (I4B) :: polylo, polyhi, nff, errcode
      real (DP) :: vert(3,4), lambda(4)
      real (DP) :: xsi(:,:), cxsi(:,:)
      logical curl
      intent (in) :: lambda, vert, polylo, polyhi, nff, curl
      intent (out) :: xsi, cxsi, errcode
      end subroutine shapefunctionv3D
end interface


interface
      subroutine shapefunctionsc3D(lambda, vert, polylo, polyhi,   &
     &  nff, grad, xsi, gxsi, errcode)
      use femtypes
      implicit none
      integer (I4B) :: polylo, polyhi, nff, errcode
      real (DP) :: vert(3,4), lambda(4)
      real (DP) :: xsi(:), gxsi(:,:)
      logical grad
      intent (in) :: lambda, vert, polylo, polyhi, nff, grad
      intent (out) :: xsi, gxsi
      intent (inout) :: errcode
      end subroutine shapefunctionsc3D
end interface


interface
      subroutine shapesc3D(l, polylo, polyhi, vec)
      Use femtypes
      implicit none
      integer (I4B) :: polylo, polyhi
      real (DP)     :: l(4)
      real (DP)     :: vec(:)
      intent (in)   :: l, polylo, polyhi
      intent (out)  :: vec
      end subroutine shapesc3D
end interface


interface
      subroutine shapev3D( l, gl, polylo, polyhi, vec )
      Use femtypes
      implicit none
      integer (I4B) :: polylo, polyhi
      real (DP)     :: l(4), gl(4)
      real (DP)     :: vec(:)
      intent (in)   :: l, polylo, polyhi
      intent (out)  :: vec
      end subroutine shapev3D
end interface


interface
      subroutine solin(ok, gilt, eleinfo, epsgl)
      use femtypes
      implicit none
      real (DP) :: epsgl
      logical :: ok, gilt, eleinfo
      intent (out) :: ok, gilt, eleinfo, epsgl
      end subroutine solin
end interface

interface
      subroutine solin_aux(ok, gilt, eleinfo, epsgl)
      use femtypes
      implicit none
      real (DP) :: epsgl
      logical :: ok, gilt, eleinfo
      intent (out) :: ok, gilt, eleinfo, epsgl
      end subroutine solin_aux
end interface

interface
      subroutine solout(gilt, eleinfo, epsgl)
      use femtypes
      implicit none
      real (DP) :: epsgl
      logical :: gilt, eleinfo
      intent (in) :: gilt, eleinfo, epsgl
      end subroutine solout
end interface


interface
      subroutine sortnodes
      use femtypes
      implicit none
      end subroutine sortnodes
end interface


interface
      subroutine SortEdges( edges , n_0 )
      use femtypes                   
      implicit none                  
      integer (I4B), pointer, optional :: n_0(:)
      integer (I4B),           pointer :: edges(:)
      intent  (inout)                  :: edges, n_0
      end subroutine SortEdges
end interface


interface
      subroutine solve_adapt3D()
      use femtypes
      implicit none
      end subroutine solve_adapt3D
end interface


interface
subroutine solve_andersone(y_n)
      use femtypes
      implicit none
      complex (DPC),      pointer :: y_n(:)
      intent (inout) :: y_n
      end subroutine solve_andersone
end interface


interface
      subroutine solve_andersonr(y_n)
      use femtypes
      implicit none
      complex (DPC),      pointer :: y_n(:)
      intent (inout) :: y_n
      end subroutine solve_andersonr
end interface


interface      
      subroutine solve_barbor(x_n)
      use femtypes
      implicit none
      complex (DPC),      pointer :: x_n(:)
      intent(inout)               :: x_n
      end subroutine solve_barbor
end interface      


interface
subroutine solve_barborv(x_n)
      use femtypes
      implicit none
      complex (DPC),      pointer :: x_n(:)
      intent(inout)               :: x_n
      end subroutine solve_barborv
end interface      


interface
      subroutine solve_diag_jacobian(x_n)
      use femtypes
      implicit none
      complex (DPC),      pointer :: x_n(:)
      intent(inout)               :: x_n
      end subroutine solve_diag_jacobian
 end interface


interface
      subroutine solve_diag_jacobianid(x_n)
      use femtypes
      implicit none
      complex (DPC),      pointer :: x_n(:)
      intent(inout)               :: x_n
      end subroutine solve_diag_jacobianid
end interface


interface
      subroutine solve_fixedpt(x_n)
      use femtypes
      implicit none
      real     (DP)               :: fp_res , fpiter_err
      complex (DPC),      pointer :: x_n(:)
      intent (inout) :: x_n
      end subroutine solve_fixedpt
end interface


interface
      subroutine solve_linear3D(epsgl, resgl, X_0, only_assembly, res_iter)
      use femtypes
      implicit none
      complex (DPC),  pointer, optional :: X_0(:), res_iter(:)
      real (DP)                         :: epsgl, resgl
      logical, optional :: only_assembly
      intent (in)                       :: X_0, only_assembly
      intent(out)                       :: epsgl, resgl, res_iter
      end subroutine solve_linear3D
end interface


interface
      subroutine solve_nonlinear3D()
      use femtypes
      implicit none
      end subroutine solve_nonlinear3D
end interface


interface
      subroutine swap23()
      use femtypes
      implicit none
      end subroutine swap23
end interface


interface
      function tetangles(vert)
      use femtypes
      implicit none
      real  (DP) :: tetangles(4), vert(3,4)
      intent(in) :: vert
      end function tetangles
end interface


interface
      pure function tetmeanratio(vert)
      USE femtypes
      implicit none
      real   (dp) :: tetmeanratio, vert(3,4)
      intent (in) :: vert
      end function tetmeanratio
end interface


interface
      pure function tetvolume( vert )
      use femtypes
      implicit none
      real   (DP) :: tetvolume, vert(3,4)
      intent (in) :: vert
      end function tetvolume
end interface


interface
      subroutine updateaux(vols,e_0,es,ev,e_len,et_nodes,surfs)
      use femtypes
      implicit none
      integer  (I4B)                    :: e_0
      integer  (I4B)                    :: vols(:), et_nodes(:)
      integer  (I4B)         , optional :: surfs(:)
      real      (DP),           pointer :: e_len(:)
      type (ARRPTRI),           pointer :: ev(:), es(:)
      intent    (in)                    :: vols, e_0, et_nodes, surfs
      intent (inout)                    :: es, ev, e_len
      end subroutine updateaux
end interface


interface
      subroutine userbc3D( face, nature, bcnum, xyzs, pval, qval )
      use femtypes
      implicit none
      integer (I4B) :: face, nature, bcnum
      real (DP) :: xyzs(3)
      complex(DPC) :: pval
      complex(DPC), optional :: qval(:)
      intent (in)  :: face, nature, bcnum, xyzs
      intent (out) :: pval, qval
      end subroutine userbc3D
end interface


interface
      subroutine userbc3D_vec( bcnum, xyzs, pval, qval )
      use femtypes
      implicit none
      integer (I4B) :: bcnum
      real (DP) :: xyzs(3)
      complex(DPC) :: pval(3)
      complex(DPC), optional :: qval(3,3)
      intent (in) :: bcnum, xyzs
      intent (out) :: pval, qval
      end subroutine userbc3D_vec
end interface


interface
      pure subroutine usermaterials3D(matname,found,xyzs,names,values,numnames, elem)
      use femtypes
      implicit none
      real (DP), optional :: xyzs(:), values(:)
      integer (I4B), optional :: numnames
      integer (I4B), optional :: elem
      character (len=*) :: matname
      character (len=*), optional :: names(:)
      logical :: found
      intent (in) :: matname, xyzs, elem
      intent (out) :: names, values, found, numnames
      end subroutine usermaterials3D
end interface


interface
      subroutine vol2aux(numf,vf)
      use femtypes
      implicit none
      integer (I4B), intent(out), optional :: numf
      integer (I4B), pointer, optional :: vf(:,:)
      end subroutine vol2aux
end interface


interface
      subroutine vtkin(ok)
      use femtypes
      implicit none
      logical ok
      intent (out) :: ok
      end subroutine vtkin
end interface

interface
      subroutine vtk_varmatplot
      use femtypes
      implicit none
      end subroutine vtk_varmatplot
end interface


interface
      subroutine vtk_scalarcellplot()
      use femtypes
      implicit none
      end subroutine vtk_scalarcellplot
end interface


interface
      subroutine vtk_scalargridplot()
      use femtypes
      implicit none
      end subroutine vtk_scalargridplot
end interface


interface
      subroutine vtk_scalarpointplot()
      use femtypes
      implicit none
      end subroutine vtk_scalarpointplot
end interface


interface
      subroutine vtk_tensorpointplot()
      use femtypes
      implicit none
      end subroutine vtk_tensorpointplot
end interface


interface
      subroutine vtk_vectorcellplot()
      use femtypes
      implicit none
      end subroutine vtk_vectorcellplot
end interface


interface
      subroutine vtk_vectorpointplot()
      use femtypes
      implicit none
      end subroutine vtk_vectorpointplot
end interface


interface
      subroutine vtk_write_mesh(unitid, ellist)
      use femtypes
      implicit none
      integer (I4B) :: unitid
      integer(I4B), optional :: ellist(:)
      intent (in) :: ellist
      end subroutine vtk_write_mesh
end interface

interface
      subroutine vtk_write_cells(unitid, cell_type, ellist)
      use femtypes
      implicit none
      integer (I4B) :: unitid, cell_type
      integer(I4B), optional :: ellist(:)
      intent (in) :: ellist, cell_type
      end subroutine vtk_write_cells
end interface


interface
      subroutine writedata(keep_olddata,filename,plottype,plottitle,xunits,yunits,xdata,ydata)
      use femtypes
      implicit none
      integer   (I4B)     :: keep_olddata
      character (len=*)   :: filename,plottype,plottitle,xunits,yunits
      real      (DP)      :: xdata(:),ydata(:)
      intent    (in)      :: filename,plottype,plottitle,xunits,yunits,xdata,ydata, keep_olddata
      end subroutine writedata
end interface


interface
      subroutine writeng(ok)
      use femtypes
      implicit none
      logical ok
      intent (out) :: ok
      end subroutine writeng
end interface


interface
      subroutine writeoutelementvalues(elements,nat)
      use femtypes
      implicit none
      real (DP)             :: elements(:)
      integer (I4B)         :: nat
      end subroutine writeoutelementvalues
end interface


interface
      subroutine writeoutnodevalues(nodes,nat)
      use femtypes
      implicit none
      real (DP)             :: nodes(:)
      integer (I4B)         :: nat
      end subroutine writeoutnodevalues
end interface


interface
      subroutine writeoutmeshmark(mark,adapt_type)
      use femtypes
      implicit none
      real (DP)             :: mark(:)
      character(len=1)      :: adapt_type
      end subroutine writeoutmeshmark
end interface


interface
      subroutine writeunv(ok)
      use femtypes
      implicit none
      logical ok
      intent (out) :: ok
      end subroutine writeunv
end interface


interface
      subroutine xyz2lam(lambda, elem, xyzt, nod, vn)
      use femtypes
      implicit none
      integer (I4B) :: elem, vn(:,:)
      real (DP) :: xyzt(3), lambda(4), nod(:,:)
      intent (in) :: xyzt, nod, elem, vn
      intent (out) :: lambda
      end subroutine xyz2lam
end interface


end module feminterface3d



