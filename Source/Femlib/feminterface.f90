module feminterface
!==========================================================================
!------------------------------------------------------------------------------
!
!    PolyDE -- a finite element simulation tool
!    Copyright (C) 2006 Institute for Micro Systems Technology,
!                       Hamburg University of Technology.
!
!    This file is part of PolyDE.
!
!    PolyDE is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by the Free Software Foundation; either version 2, or (at your
!    option) any later version.
!
!    PolyDE is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!    MA 02110-1301 USA.
!
!    $Revision: 1.159 $
!    $Date: 2015/09/11 15:15:24 $
!    $Author: juryanatzki $
!
!------------------------------------------------------------------------------
interface a_norm
      function a_norm_DD(x,diag,lower,upper,lia,lja)
      use femtypes
      implicit none
      real(DP) :: diag(:), lower(:), upper(:)
      real(DP) :: x(:)
      real(DP) :: a_norm_DD
      integer(I4B) :: lia(:), lja(:)
      intent(in) :: x, diag, lower, upper, lia, lja
      end function a_norm_DD

      function a_norm_SS(x,diag,lower,upper,lia,lja)
      use femtypes
      implicit none
      real(SP) :: diag(:), lower(:), upper(:)
      real(SP) :: x(:)
      real(SP) :: a_norm_SS
      integer(I4B) :: lia(:), lja(:)
      intent(in) :: x, diag, lower, upper, lia, lja
      end function a_norm_SS
end interface a_norm


interface
      subroutine abline(x1,y1,x2,y2,x3,y3,x4,y4,abst,len1,len2)
      use femtypes
      implicit none
      real (DP) x1,y1,x2,y2,x3,y3,x4,y4,abst,len1,len2
      intent (in) :: x1, y1, x2, y2, x3, y3, x4, y4
      intent (out) :: abst, len1, len2
      end subroutine abline
end interface


interface
      subroutine acadini
      use femtypes
      implicit none
      end subroutine acadini
end interface


interface
      subroutine adaptation(epsgl)
      use femtypes
      implicit none
      real (DP), intent (out) :: epsgl
      end subroutine adaptation
end interface


interface
      subroutine amux ( n, x, y, a, ja, ia )
      use femtypes
      implicit none
      integer (I4B) :: n, ia(*), ja(*)
      complex (DPC) ::  a(*), x(*), y(n)
      intent (in) :: n, a, ia, ja, x
      intent (out) :: y
      end subroutine amux
end interface


interface
      subroutine angles(i,e,xn,yn,w1,w2,w3)
      use femtypes
      implicit none
      integer (I4B) i, e(:,:)
      real (DP) xn(:), yn(:), w1, w2, w3
      intent (in) :: i, e, xn, yn
      intent (out) w1, w2, w3
      end subroutine angles
end interface


interface
      subroutine arc(accur,fakt1,scalfk,laytxt,layrb,layanz,xpoint,     &
     &  ypoint,zki,lzrb,zpz,length,anzzwg,anzknt,ok,txtlen)
      use femtypes
      implicit none
      integer (I4B) :: anzzwg, anzknt, txtlen, layanz
      integer (I4B), pointer :: lzrb(:), zki(:,:), zpz(:), layrb(:,:)
      real (DP) :: accur, fakt1, scalfk
      real (DP), pointer :: xpoint(:), ypoint(:), length(:)
      logical :: ok
      character (len=txtlen), pointer :: laytxt(:)
      intent (in) :: accur, fakt1, scalfk
      intent (out) :: xpoint, ypoint, zki, lzrb, zpz, length, ok
      intent (inout) :: laytxt, layrb, layanz, anzzwg, anzknt
      end subroutine arc
end interface


interface
      subroutine assembly(lower,upper,diag,acsr,rhs,ia,ja,jacobi,csr,   &
     &                      matvar,symmetric)
      use femtypes
      implicit none
      integer (I4B), pointer :: ia(:), ja(:)
      complex (DPC), pointer :: lower(:), upper(:), diag(:), rhs(:), acsr(:)
      logical jacobi, csr, matvar, symmetric
      intent (in) ::  jacobi, csr, matvar, symmetric
      end subroutine assembly
end interface


interface
      subroutine assemblyandsolve(epsgl,resgl,adaptstep)
      use femtypes
      implicit none
      integer (I4B), optional :: adaptstep
      real (DP)   :: epsgl, resgl
      intent (in) :: adaptstep
      intent (out) :: epsgl, resgl
      end subroutine assemblyandsolve
end interface


interface
      subroutine assemblyandsolvetrans(epsgl,resgl)
      use femtypes
      implicit none
      real (DP)   :: epsgl, resgl
      intent (out) :: epsgl, resgl
      end subroutine assemblyandsolvetrans
end interface


interface
      subroutine assemblytrans(lower,upper,diag,acsrs,acsrd,acsrm,rhs,ia,ja,jacobi,csr,   &
     &                      matvar,symmetric)
      use femtypes
      implicit none
      integer (I4B), pointer :: ia(:), ja(:)
      complex (DPC), pointer :: lower(:), upper(:), diag(:), acsrs(:), acsrd(:), acsrm(:), rhs(:)
      logical jacobi, csr, matvar, symmetric
      intent (in) ::  jacobi, csr, matvar, symmetric
      end subroutine assemblytrans
end interface


interface atxb
      subroutine atxb_DPC(lower,upper,diag,b,x,ia,ja,n)
      use femtypes
      implicit none
      complex (DPC) lower(:),upper(:),diag(:)
      complex (DPC) x(:),b(:)
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
      end subroutine atxb_DPC

      subroutine atxb_DP(lower,upper,diag,b,x,ia,ja,n)
      use femtypes
      implicit none
      real (DP) lower(:),upper(:),diag(:)
      real (DP) x(:),b(:)
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
      end subroutine atxb_DP
end interface


interface
      subroutine ausgab(ele,xk,yk,nb,ez,kz,kzi,zrb,krb,geb,konfom,zyl)
      use femtypes
      implicit none
      integer (I4B) ele(:,:),nb(:,:),ez,kz,kzi(:),zrb(:,:),krb(:,:), geb(:)
      real (DP) xk(:),yk(:)
      logical konfom,zyl
      intent (in) :: nb, ez, kz, zrb, krb, geb, konfom, zyl
      intent (inout) :: ele, xk, yk, kzi
      end subroutine ausgab
end interface


interface axb
      subroutine axb_DPC(lower,upper,diag,b,x,ia,ja,n)
      use femtypes
      implicit none
      complex (DPC) lower(:),upper(:),diag(:)
      complex (DPC) x(:),b(:)
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
      end subroutine axb_DPC

      subroutine axb_DP(lower,upper,diag,b,x,ia,ja,n)
      use femtypes
      implicit none
      real (DP) lower(:),upper(:),diag(:)
      real (DP) x(:),b(:)
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
      end subroutine axb_DP
end interface


interface
      subroutine bandw(e,n,bd,liste,sip)
      use femtypes
      implicit none
      integer (I4B) e(:,:),n,bd
      integer (I4B) liste(:),sip
      intent (in) :: e, n, liste, sip
      intent (out) :: bd
      end subroutine bandw
end interface


interface
      subroutine basin(ok,jd,n,p,e,xn,yn,kzi,geb,en)
      use femtypes
      implicit none
      integer (I4B) n, p 
      integer (I4B), pointer :: e(:,:), kzi(:), geb(:), en(:,:)
      real (DP) :: jd
      real (DP), pointer :: xn(:), yn(:)
      logical ok
      intent (out) ok, n, p, jd
      end subroutine basin
end interface


interface
      subroutine basout(ok,n,p,e,xn,yn,kzi,geb,en)
      use femtypes
      implicit none
      integer (I4B) e(:,:),n,p,kzi(:),geb(:),en(:,:)
      real (DP) xn(:),yn(:)
      logical ok
      intent (in) :: n, p, e, xn, yn, kzi, geb, en
      intent (out) :: ok
      end subroutine basout
end interface


interface
      subroutine bfarb(a,zcolmin,zcolmax,x,y,e,ie,h)
      use femtypes
      implicit none
      real (DP) zcolmin,zcolmax,x(:),y(:),h
      complex (DPC) a(:)
      integer (I4B) e(:,:),ie
      intent (in):: a,zcolmin,zcolmax,x,y,e,ie,h
      end subroutine bfarb
end interface


interface
      subroutine bfarbi(a,zcolmin,zcolmax,x,y,e,ie,h,p)
      use femtypes
      implicit none
      real (DP) zcolmin,zcolmax,x(:),y(:),h
      complex (DPC) a(:)
      integer (I4B) e(:,:),ie,p
      intent (in):: a,zcolmin,zcolmax,x,y,e,ie,h,p
      end subroutine bfarbi
end interface


interface
      subroutine bh(b,myr,dnyb2)
      use femtypes
      implicit none
      real (DP) b,myr,dnyb2
      intent (in) :: b
      intent (out) :: myr, dnyb2
      end subroutine bh
end interface


interface
      subroutine blocks
      use femtypes
      implicit none
      end subroutine blocks
end interface


interface
      subroutine bmaxr(x,e,n,xn,yn)
      use femtypes
      implicit none
      real (DP) xn(:), yn(:)
      complex (DPC) x(:)
      integer (I4B) e(:,:), n
      intent (in) :: x, e, n, xn, yn
      end subroutine bmaxr
end interface


interface bsor
      subroutine bsor_d(upper,ia,ja,b,x,n,om1)
      use femtypes
      implicit none
      real (DP) upper(:),b(:),x(:),om1
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: upper, ia, ja, n, om1
      intent (out) :: x
      intent (inout) :: b
      end subroutine bsor_d

      subroutine bsor_c(upper,ia,ja,b,x,n,om1)
      use femtypes
      implicit none
      complex (SPC) upper(:),b(:),x(:)
      real (SP) om1
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: upper, ia, ja, n, om1
      intent (out) :: x
      intent (inout) :: b
      end subroutine bsor_c

      subroutine bsor_z(upper,ia,ja,b,x,n,om1)
      use femtypes
      implicit none
      complex (DPC) upper(:),b(:),x(:)
      real (DP) om1
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: upper, ia, ja, n, om1
      intent (out) :: x
      intent (inout) :: b
      end subroutine bsor_z

      subroutine bsor_s(upper,ia,ja,b,x,n,om1)
      use femtypes
      implicit none
      real (SP) upper(:),b(:),x(:),om1
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: upper, ia, ja, n, om1
      intent (out) :: x
      intent (inout) :: b
      end subroutine bsor_s
end interface


interface
      subroutine cal2juld(jd,year,mon,day,hour,minute,second)
      use femtypes
      implicit none
      real (DP) jd
      integer (I4B) year, mon, day, hour, minute, second
      intent (in) year, mon, day, hour, minute, second
      intent (out) jd
      end subroutine cal2juld
end interface


interface calctensor
      pure subroutine calctensor_c(var11, var22, angle, invert, ten)
      use femtypes
      implicit none
      real (DP) :: angle
      complex (DPC) :: var11, var22, ten(2,2)
      logical :: invert
      intent (in) :: angle, var11, var22, invert
      intent (out) :: ten
      end subroutine calctensor_c

      pure subroutine calctensor_d(var11, var22, angle, invert, ten)
      use femtypes
      implicit none
      real (DP) :: angle
      real (DP) :: var11, var22, ten(2,2)
      logical :: invert
      intent (in) :: angle, var11, var22, invert
      intent (out) :: ten
      end subroutine calctensor_d
end interface


interface
      subroutine calctimestep_newmark(d2xdt2_at_tminusdt,d2xdt2_at_tplusdt,timestep,nexttimestep, &
      &          currenttime,alpha_newmark,tranltollo,tranltolup,timeglobalerror)
      use femtypes, only: DP, DPC
      implicit none
      complex (DPC), pointer :: d2xdt2_at_tminusdt(:), d2xdt2_at_tplusdt(:)
      real (DP) :: timestep, nexttimestep, currenttime, alpha_newmark
      real (DP) :: tranltollo, tranltolup, timeglobalerror
      intent (in) :: timestep, currenttime, alpha_newmark, tranltollo, tranltolup
      intent (out) :: nexttimestep, timeglobalerror
      end subroutine calctimestep_newmark
end interface


interface
      subroutine calctimestep_runge(k1,k2,timestep,nexttimestep,timestepreject)
      use femtypes
      implicit none
      complex (DPC), pointer :: k1(:), k2(:)
      real (DP)   :: timestep, nexttimestep
      logical :: timestepreject
      intent (in) :: timestep
      intent (out) :: nexttimestep
      intent (inout) :: timestepreject
      end subroutine calctimestep_runge
end interface


interface
      subroutine centerofgravity(x,y,n,xctr,yctr)
      use femtypes
      implicit none
      real (DP) x(:), y(:), xctr, yctr
      integer (I4B) n
      intent (in) :: x, y, n
      intent (out) :: xctr, yctr
      end subroutine centerofgravity
end interface


interface
      subroutine cgauss(ca,cb,n,cx)
      use femtypes
      implicit none
      integer (I4B) n
      complex (DPC) ca(:,:),cb(:),cx(:)
      intent (in) :: ca,cb,n
      intent (out) :: cx
      end subroutine cgauss
end interface


interface
      subroutine classes
      use femtypes
      implicit none
      end subroutine classes
end interface


interface
      subroutine circ(accur,fakt1,scalfk,laytxt,layrb,layanz,xpoint,    &
     &  ypoint,zki,lzrb,zpz,length,anzzwg,anzknt,ok,txtlen)
      use femtypes
      implicit none
      integer (I4B) :: anzzwg, anzknt, txtlen, layanz
      integer (I4B), pointer :: lzrb(:), zki(:,:), zpz(:), layrb(:,:)
      real (DP) :: accur, fakt1, scalfk
      real (DP), pointer :: xpoint(:), ypoint(:), length(:)
      logical :: ok
      character (len=txtlen), pointer :: laytxt(:)
      intent (in) :: accur, fakt1, scalfk
      intent (out) :: xpoint, ypoint, zki, lzrb, zpz, length, ok
      intent (inout) :: laytxt, layrb, layanz, anzzwg, anzknt
      end subroutine circ
end interface


interface
      subroutine circle(x, y, w1, w2, r1, r2)
      use femtypes
      implicit none
      real (SP) x, y, w1, w2, r1, r2
      intent (in)::  x, y, w1, w2, r1, r2
      end subroutine circle
end interface


interface coloring
      subroutine coloring(indx,offset)
      use femtypes
      integer (I4B), allocatable :: indx(:), offset(:)
      intent (out) :: indx, offset
      end subroutine coloring
end interface coloring


interface
      subroutine concat(plist,plist2,plp,plp2,merk,rand,                &
     &  lauf,zaehl,zahl,pzr)
      use femtypes
      implicit none
      integer (I4B) plist(:),plist2(:),plp(:),plp2(:)
      integer (I4B) merk,rand,lauf,zaehl,zahl,pzr
      intent (in) :: merk, rand, lauf, zaehl
      intent (out) :: plist2, plp2, pzr
      intent (inout) :: plist, plp, zahl
      end subroutine concat
end interface


interface
      subroutine convertcolorstring(red_f, green_f, blue_f)
      use femtypes
      implicit none
      real (SP) :: red_f, green_f, blue_f
      intent (out) :: red_f, green_f, blue_f
      end subroutine convertcolorstring
end interface


interface
      subroutine copyr(text1,text2)
      use femtypes
      implicit none
      character (len=*) text1,text2
      intent (in) :: text1, text2
      end subroutine copyr
end interface


interface cross_product
      pure function cross_product_i(a,b)
      use femtypes
      implicit none
      integer (I4B) :: a(3), b(3), cross_product_i(3)
      intent (in) :: a, b
      end function cross_product_i

      pure function cross_product_r(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3), b(3), cross_product_r(3)
      intent (in) :: a, b
      end function cross_product_r

      pure function cross_product_dpc(a,b)
      use femtypes
      implicit none
      complex (DPC) :: a(3), b(3), cross_product_dpc(3)
      intent (in) :: a, b
      end function cross_product_dpc

      pure function cross_product_rdpc(a,b)
      use femtypes
      implicit none
      complex (DPC) :: b(3), cross_product_rdpc(3)
      real (DP)     :: a(3)
      intent (in) :: a, b
      end function cross_product_rdpc

      pure function cross_product_dpcr(a,b)
      use femtypes
      implicit none
      complex (DPC) :: a(3), cross_product_dpcr(3)
      real (DP)     :: b(3)
      intent (in) :: a, b
      end function cross_product_dpcr

end interface


interface csrcsc
      subroutine csrcsc_z( n, a, ja, ia, ao, jao, iao )
      use femtypes
      implicit none
      integer (I4B) :: n, ia(n+1), iao(n+1), ja(*), jao(*)
      complex (DPC) ::  a(*), ao(*)
      intent (in) :: n, a, ia, ja
      intent (out) :: ao, iao, jao
      end subroutine csrcsc_z

      subroutine csrcsc_d( n, a, ja, ia, ao, jao, iao )
      use femtypes
      implicit none
      integer (I4B) :: n, ia(n+1), iao(n+1), ja(*), jao(*)
      real (DP) ::  a(*), ao(*)
      intent (in) :: n, a, ia, ja
      intent (out) :: ao, iao, jao
      end subroutine csrcsc_d
end interface


interface csrmmout
    subroutine csrmmout(ia,ja,csr,adaptstep)
    use femtypes
    implicit none
    integer (I4B), pointer :: ia(:) ,ja(:)
    complex (DPC), pointer :: csr(:)
    integer (I4B), optional :: adaptstep
    intent (in) :: ia, ja, csr, adaptstep
    end subroutine csrmmout
end interface


interface
      subroutine currentflow(xa,ya,xe,ye,current)
      use femtypes
      implicit none
      real (DP) :: xa, ya, xe, ye
      real (DP) :: current
      intent (in) :: xa, ya, xe, ye
      intent (out) :: current
      end subroutine currentflow
end interface


interface
      subroutine decroach(encrlimit, meshsiz2)
      use femtypes
      implicit none
      real (DP) encrlimit, meshsiz2(:)
      intent (in) :: encrlimit, meshsiz2
      end subroutine decroach
end interface


interface
      subroutine decroachn(encrlimit,x3,y3,splitted)
      use femtypes
      implicit none
      real (DP) encrlimit, x3, y3
      logical splitted
      intent (in) :: encrlimit, x3, y3
      intent (out) :: splitted
      end subroutine decroachn
end interface


interface
      subroutine degree(root, xadj, adjncy, mask, deg, ccsize)
      use femtypes
      implicit none
      integer (I4B) root, xadj(:), adjncy(:), mask(:), deg(:), ccsize
      intent (in) :: root, adjncy, mask
      intent (out) :: deg, ccsize
      intent (inout) :: xadj
      end subroutine degree
end interface


interface
      subroutine delaunay_refinement(meshsize)
      use femtypes
      implicit none
      real (DP), pointer:: meshsize(:)
      end subroutine delaunay_refinement
end interface


interface destroyarrptr
      function destroyarrptr_i(arrptr)
      use femtypes
      implicit none
      type (ARRPTRI), pointer :: arrptr(:)
      logical :: destroyarrptr_i
      end function destroyarrptr_i

      function destroyarrptr_i2(arrptr)
      use femtypes
      implicit none
      type(ARRPTRI), pointer :: arrptr(:,:)
      logical :: destroyarrptr_i2
      end function destroyarrptr_i2

      function destroyarrptr_r(arrptr)
      use femtypes
      implicit none
      type (ARRPTRR), pointer :: arrptr(:)
      logical :: destroyarrptr_r
      end function destroyarrptr_r

      function destroyarrptr_dpc(arrptr)
      use femtypes
      implicit none
      type (ARRPTRDPC), pointer :: arrptr(:)
      logical :: destroyarrptr_dpc
      end function destroyarrptr_dpc
end interface


interface
      pure subroutine dgshape(l, nu, xl, yl, polylo, polyhi, nff, vec, errcode)
      use femtypes
      implicit none
      integer (I4B) elem, polylo, polyhi, nff, errcode
      real (DP) xl(:), yl(:), l(:)
      complex (DPC) vec(:), nu(:,:)
      intent (in) :: l, nu, xl, yl, polylo, polyhi, nff
      intent (out) :: vec, errcode
      end subroutine dgshape
end interface


interface
      subroutine dichk(n,p,e,en,xn,yn,ok,geb,kzi)
      use femtypes
      implicit none
      integer (I4B) n,p,e(:,:),en(:,:)
      integer (I4B), optional :: geb(:), kzi(:)
      real (DP) xn(:),yn(:)
      logical ok
      intent (in) :: n, p, e, xn, yn, geb, kzi
      intent (out) :: ok
      intent (inout) :: en
      end subroutine dichk
end interface


interface
      subroutine divtri(xl,yl,zl,zval,first,all,done,npkt,xpoly,ypoly,zpoly)
      use femtypes
      implicit none
      real (DP) xl(3), yl(3), zl(3), zval
      real (SP) xpoly(:), ypoly(:), zpoly(:)
      integer (I4B) npkt
      logical done, first, all
      intent (in):: xl, yl, zl, zval, all
      intent (inout):: first
      intent (out) done, npkt, xpoly, ypoly, zpoly
      end subroutine divtri
end interface


interface
      subroutine dptransformpoly(xl,yl,zl,a,n,xo,yo,zo)
      use femtypes
      implicit none
      real (DP) a(4,4), xl(:), yl(:), zl(:), xo(:), yo(:), zo(:)
      integer (I4B) n
      intent (in) xl, yl, zl, a, n
      intent (out) xo, yo, zo
      end subroutine dptransformpoly
end interface


interface
      subroutine drtri(i,e,x,y,z,zcolor,ncolor,palette,proj,xo,yo,zo,transf)
      use femtypes
      implicit none
      real (DP) x(:),y(:),z(:),zcolor(:)
      real (DP), optional ::  xo(:),yo(:),zo(:),transf(4,4)
      integer (I4B) i,e(:,:),ncolor,palette(:)
      logical proj
      intent (in):: i,e,x,y,z,xo,yo,zo,zcolor,ncolor,palette,transf,proj
      end subroutine drtri
end interface


interface
      subroutine earclipper(head,numnodes,region)
      use femtypes
      use triangtypes
      implicit none
      integer (I4B) region, numnodes
      type (node), pointer :: head
      intent (in) :: region
      intent (inout) :: numnodes
      end subroutine earclipper
end interface


interface
      subroutine eartest(minode,head,nbds,actbd,xn,yn,numnodes,ear)
      use femtypes
      use triangtypes
      implicit none
      integer (I4B) nbds, actbd, numnodes(:)
      real (DP) xn(:), yn(:)
      logical ear
      type (node), pointer :: minode
      type (arrpt), pointer :: head(:)
      intent (in) :: nbds, actbd, xn, yn, numnodes
      intent (out) :: ear
      end subroutine eartest
end interface


interface
      subroutine elementmatrix(elem,polyorder,jacobi,full,matvar,a,b,nff,nffsum,errcode)
      use femtypes
      implicit none
      integer (I4B) :: elem
      integer (I4B) :: errcode
      integer (I4B) :: polyorder(:)
      integer (I4B), dimension(:), allocatable :: nff, nffsum
      complex(DPC), pointer :: a(:,:), b(:)
      logical :: jacobi, full, matvar
      intent (in) :: elem, polyorder, jacobi, full, matvar
    !  intent (out) :: nff, nffsum  ! org
      intent (inout) :: nff, nffsum ! test
      end subroutine elementmatrix
end interface


interface
      subroutine elementmatrixtrans(elem,polyorder,jacobi,full,matvar,s,d,m,b,nff,nffsum,errcode)
      use femtypes
      implicit none
      integer (I4B) :: elem
      integer (I4B) polyorder(:), nff(:), nffsum(:), errcode
      complex (DPC), pointer :: s(:,:), d(:,:), m(:,:), b(:)
      logical :: jacobi, full, matvar
      intent (in) :: elem, polyorder, jacobi, full, matvar
      intent (out) :: nff,nffsum
      end subroutine elementmatrixtrans
end interface


interface
      subroutine elemnd(xn,yn,e,n,xt,yt,ielem,ok)
      use femtypes
      implicit none
      real (DP) xn(:), yn(:), xt, yt
      integer (I4B) e(:,:), n, ielem
      logical ok
      intent (in) :: xn, yn, e, n, xt, yt
      intent (out) :: ielem, ok
      end subroutine elemnd
end interface


interface
      subroutine elemnt(xn,yn,e,n,xt,yt,ielem,en,ok)
      use femtypes
      implicit none
      real (DP) xn(:), yn(:), xt, yt
      integer (I4B) e(:,:), n, ielem, en(:,:)
      logical ok
      intent (in) :: xn, yn, e, n, xt, yt, en
      intent (out) :: ielem, ok
      end subroutine elemnt
end interface


interface        
      subroutine elemtx(i,x,a,b,jacobi,xn,yn,e,omega,matzif,geb)
      use femtypes
      implicit none
      integer (I4B) i,e(:,:),matzif(:),geb(:)
      real (DP) xn(:),yn(:),omega
      complex (DPC) a(3,3),b(3),x(:)
      logical jacobi
      intent (in) :: i, x, jacobi, xn, yn, e, omega, matzif, geb
      intent (out) :: a,b
      end subroutine elemtx
end interface


interface
      subroutine elmgeb(gbz,bzi,bzip,bzil,lzrb,layrb,matname)
      use femtypes
      implicit none
      integer (I4B) :: gbz, bzi(:), bzip(:), bzil(:), lzrb(:), layrb(:,:)
      character (len=*) :: matname(:)
      intent (in) :: lzrb, layrb, matname
      intent (inout) :: gbz, bzi, bzip, bzil
      end subroutine elmgeb
end interface


interface
      subroutine encrtest(i,j,x3,y3,encrlimit,e,en,geb,kzi,xn,yn,encroached)
      use femtypes
      implicit none
      integer (I4B) i, j, e(:,:),  en(:,:), geb(:), kzi(:)
      real (DP) encrlimit, xn(:), yn(:), x3, y3
      logical encroached
      intent (in) :: i, j, x3, y3, encrlimit, e, en, geb, kzi, xn, yn
      intent (out) :: encroached
      end subroutine encrtest
end interface


interface
      subroutine energyflow(xa,ya,xe,ye,flux)
      use femtypes
      implicit none
      real (DP) :: xa, ya, xe, ye
      real (DP) :: flux
      intent (in) :: xa, ya, xe, ye
      intent (out) :: flux
      end subroutine energyflow
end interface


interface
      subroutine entity(accur,fakt1,scalfk,xpoint,ypoint,zki,lzrb,zpz,  &
     &  xgeb,ygeb,txtgeb,length,anzzwg,anzknt,anzgeb,laytxt,layrb,      &
     &  lplayt,layanz,snapof,txtlen)
      use femtypes
      implicit none
      integer (I4B) :: anzzwg, anzknt, anzgeb, snapof, layanz, txtlen
      integer (I4B), pointer :: lzrb(:), zki(:,:), zpz(:), lplayt(:), layrb(:,:)
      real (DP) :: accur, fakt1, scalfk
      real (DP), pointer :: xgeb(:), ygeb(:), xpoint(:), ypoint(:), length(:)
      character (len=txtlen), pointer :: txtgeb(:), laytxt(:)
      intent (in) :: accur, fakt1, scalfk, txtlen
      intent (out) :: xpoint, ypoint, zki, lzrb, zpz
      intent (out) :: length, anzzwg, anzknt, anzgeb, lplayt, snapof
      intent (inout) :: laytxt, layrb, layanz, txtgeb, xgeb, ygeb
      end subroutine entity
end interface


interface
      subroutine error(sectin,level,errnum,errtxt)
      use femtypes
      implicit none
      integer (I4B) level,errnum
      character (len=*) errtxt
      character (len=*) sectin
      intent (in) :: sectin, level, errnum, errtxt
      end subroutine error
end interface


interface
      subroutine extend(anzzwg,xpoint,ypoint,zki,                       &
     &                  xmin,xmax,ymin,ymax)
      use femtypes
      implicit none
      real (DP) xpoint(:),ypoint(:),xmin,xmax,ymin,ymax
      integer (I4B) anzzwg,zki(:,:)
      intent (in) anzzwg,xpoint,ypoint,zki 
      intent (out) xmin,xmax,ymin,ymax
      end subroutine extend
end interface


interface
      subroutine farbe(z,x,y,e,ie,palette,ncolor,zcolmin,zcolmax)
      use femtypes
      implicit none
      real (DP) x(:), y(:), z(:), zcolmin, zcolmax
      integer (I4B) e(:,:), ie, ncolor, palette(:)
      intent (in):: z, x, y, e, ie, palette, ncolor, zcolmin, zcolmax
      end subroutine farbe
end interface


interface f_norm
      function f_norm_DD(diag,lower,upper,lia,lja)
      use femtypes
      implicit none
      real(DP) :: diag(:), lower(:), upper(:)
      real(DP) :: f_norm_DD
      integer(I4B) :: lia(:), lja(:)
      intent(in) :: diag, lower, upper, lia, lja
      end function f_norm_DD

      function f_norm_SS(diag,lower,upper,lia,lja)
      use femtypes
      implicit none
      real(SP) :: diag(:), lower(:), upper(:)
      real(SP) :: f_norm_SS
      integer(I4B) :: lia(:), lja(:)
      intent(in) :: diag, lower, upper, lia, lja
      end function f_norm_SS
end interface f_norm


interface
      subroutine farbpalette(palette,ncolor,color)
      use femtypes
      implicit none
      integer (I4B) ncolor
      integer (I4B) palette(:)
      character (len=*) :: color
      intent (in) ::ncolor, color
      intent (out):: palette
      end subroutine farbpalette
end interface


interface
      pure subroutine fetchmatparameters(list,matindex,xyzs,elem)
      use femtypes
      implicit none
      integer (I4B) matindex
      integer (I4B), optional :: elem
      real(DP) list(:), xyzs(:)
      intent (in) :: matindex, xyzs, elem
      intent (out) :: list
      end subroutine fetchmatparameters
end interface


interface
      subroutine ffein(ellist,zahle,zahlk,ende,merror,maxelemt,ext,int)
      use femtypes
      implicit none
      integer (I4B) :: ellist(:), zahle, zahlk, maxelemt
      real (DP) :: merror
      logical :: ende, ext, int
      intent (in) :: maxelemt, ext, int
      intent (out) :: ende, merror 
      intent (inout) :: ellist, zahle, zahlk
      end subroutine ffein
end interface



interface
      subroutine ffeinhp(ellist,zahle,zahlk,ende,merror,maxelemt,ext,int)
      use femtypes
      implicit none
      integer (I4B) ellist(:), zahle, zahlk, maxelemt
      real (DP) merror
      logical :: ende, ext, int
      intent (in) :: maxelemt, ext, int
      intent (out) :: ende, merror
      intent (inout) :: ellist, zahle, zahlk
      end subroutine ffeinhp
end interface


interface
      subroutine ffeinhp_heuveline(refinelist,resold,ellist,zahle,zahlk,ende,error,maxelemt,ext,int,res)
      use femtypes
      implicit none
      integer (I4B) ellist(:), refinelist(:), zahle, zahlk, maxelemt
      real (DP) error, resold(:), res(:)
      logical :: ende, ext, int
      intent (in) :: maxelemt, ext, int, resold
      intent (out) :: ende, error, res
      intent (inout) :: ellist, zahle, zahlk, refinelist
      end subroutine ffeinhp_heuveline
end interface


interface
      subroutine ffeinhp_KPphas(ellist,zahle,zahlk,ende,merror,maxelemt,ext,int)
      use femtypes
      implicit none
      integer (I4B) ellist(:), zahle, zahlk, maxelemt
      real (DP) merror
      logical :: ende, ext, int
      intent (in) :: maxelemt, ext, int
      intent (out) :: ende, merror 
      intent (inout) :: ellist, zahle, zahlk
      end subroutine ffeinhp_KPphas
end interface


interface
      subroutine ffeinhp_melenk(prederror,refinelist,ellist,zahle,zahlk,ende,error,maxelemt,ext,int,res)
      use femtypes
      implicit none
      integer (I4B) ellist(:), refinelist(:), zahle, zahlk, maxelemt
      real (DP) error, prederror(:), res(:)
      logical :: ende, ext, int
      intent (in) :: maxelemt, ext, int, prederror
      intent (out) :: ende, error, res
      intent (inout) :: ellist, zahle, zahlk, refinelist
      end subroutine ffeinhp_melenk
end interface


interface
      subroutine ffeinhp_top5(ellist,zahle,zahlk,ende,merror,maxelemt,ext,int)
      use femtypes
      implicit none
      integer (I4B) ellist(:), zahle, zahlk, maxelemt
      real (DP) merror
      logical :: ende, ext, int
      intent (in) :: maxelemt, ext, int
      intent (out) :: ende, merror 
      intent (inout) :: ellist, zahle, zahlk
      end subroutine ffeinhp_top5
end interface


interface
      subroutine field(elem,lambda,typ,z,soln,u,alphau,gammau,gradu,betagradu,gammagradu,nugradu,dgradu,f,g)
      use femtypes
      implicit none
      integer (I4B) elem
      real (DP) lambda(3)
      logical typ(5)
      complex (DPC) z(:,:)      !  second index for nature
      complex (DPC), optional :: soln(:)
      complex (DPC), optional :: u(:), alphau(:), betagradu(:), gammagradu(:),  dgradu(:), f(:)
      complex (DPC), optional :: gammau(:,:), gradu(:,:), nugradu(:,:), g(:,:) 
      intent (in) ::  elem, lambda, typ, soln
      intent (out) ::  z, u, alphau, gammau, gradu, betagradu, gammagradu, nugradu, dgradu, f, g
      end subroutine field
end interface


interface fieldquantity
      subroutine fieldquantity_scalar(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
      use femtypes
      implicit none
      real (DP) :: lambda(3), phi
      complex (DPC) :: zs
      character (len=*) :: fieldtype, unit
      character (len=*) :: descriptor
      integer (I4B) :: elem
      logical :: ok
      intent (in) :: elem, fieldtype, lambda, phi
      intent (out) :: descriptor, unit, zs, ok
      end subroutine fieldquantity_scalar

      subroutine fieldquantity_vector(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
      use femtypes
      implicit none
      real (DP) :: lambda(3), phi
      complex (DPC) :: zs(:)
      character (len=*) :: fieldtype, unit
      character (len=*) :: descriptor
      integer (I4B) :: elem
      logical :: ok
      intent (in) :: elem, fieldtype, lambda, phi
      intent (out) :: descriptor, unit, zs, ok
      end subroutine fieldquantity_vector

      subroutine fieldquantity_tensor(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
      use femtypes
      implicit none
      real (DP) :: lambda(3), phi
      complex (DPC) :: zs(:,:)
      character (len=*) :: fieldtype, unit
      character (len=*) :: descriptor
      integer (I4B) :: elem
      logical :: ok
      intent (in) :: elem, fieldtype, lambda, phi
      intent (out) :: descriptor, unit, zs, ok
      end subroutine fieldquantity_tensor
end interface


interface
      pure subroutine fieldu(elem, lambda, u, gradu)
      use femtypes
      implicit none
      integer (I4B) elem
      real (DP) lambda(3)
      complex (DPC) :: u(:)
      complex (DPC), optional :: gradu(:,:)
      intent (in) ::  elem, lambda
      intent (out) ::  u, gradu
      end subroutine fieldu
end interface


interface
      pure real (DP) function flaech(x1,y1,x2,y2,x3,y3)
      use femtypes
      implicit none
      real (DP) x1,y1,x2,y2,x3,y3
      intent (in) :: x1,y1,x2,y2,x3,y3
      end function flaech
end interface


interface
      subroutine fnroot(root, xadj, adjncy, mask, nlvl, ls, offset)
      use femtypes
      implicit none
      integer (I4B) root, xadj(:), adjncy(:), mask(:), nlvl
      integer (I4B) ls(:), offset
      intent (in) :: xadj, adjncy, offset
      intent (out) :: nlvl, ls
      intent (inout) :: root, mask
      end subroutine fnroot
end interface


interface
      subroutine four1(data,nn,isign)
      use femtypes
      integer (I4B) nn, isign
      real (DP) data(*)
      end subroutine four1
end interface


interface
      subroutine fpotential(a,ie,x,y,e,p,nc,phi,amax,amin)
      use femtypes
      implicit none
      real (DP) x(:), y(:), phi, amax, amin
      complex (DPC) a(:)
      integer (I4B) e(:,:), ie, p, nc
      intent (in):: a, ie, x, y, e, p, nc, phi, amax, amin
      end subroutine fpotential
end interface


interface
      subroutine fpotentialp(nc,nl,phi,amax,amin,colorfill,hlin,        &
     &                       r,g,b,fieldtype)
      use femtypes
      implicit none
      real (SP) :: r, g, b
      real (DP) :: phi, amax, amin
      integer (I4B) :: nc, nl
      logical :: colorfill, hlin
      character (len=*) :: fieldtype
      intent (in):: nc, nl, phi, amax, amin, colorfill, hlin
      intent (in):: r, g, b, fieldtype
      end subroutine fpotentialp
end interface


interface fsor
      subroutine fsor_c(lower,ia,ja,b,x,n,om1)
      use femtypes
      implicit none
      complex (SPC) lower(:),x(:),b(:)
      real (SP) om1
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: lower, ia, ja, b, n, om1
      intent (out) :: x
      end subroutine fsor_c

      subroutine fsor_z(lower,ia,ja,b,x,n,om1)
      use femtypes
      implicit none
      complex (DPC) lower(:),x(:),b(:)
      real (DP) om1
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: lower, ia, ja, b, n, om1
      intent (out) :: x
      end subroutine fsor_z

      subroutine fsor_d(lower,ia,ja,b,x,n,om1)
      use femtypes
      implicit none
      real (DP) lower(:),x(:),b(:),om1
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: lower, ia, ja, b, n, om1
      intent (out) :: x
      end subroutine fsor_d

      subroutine fsor_s(lower,ia,ja,b,x,n,om1)
      use femtypes
      implicit none
      real (SP) lower(:),x(:),b(:),om1
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: lower, ia, ja, b, n, om1
      intent (out) :: x
      end subroutine fsor_s
end interface


interface
      subroutine fsortn(felst,liste,n)
      use femtypes
      implicit none
      integer (I4B) n,liste(:)
      real (DP) felst(:)
      intent (in) :: felst,n
      intent (out):: liste
      end subroutine fsortn
end interface


interface
      subroutine ftlay(string,layanz,laytxt,layrb,layer,txtlen)
      use femtypes
      implicit none
      integer (I4B) :: layanz, layer, txtlen
      integer (I4B), pointer :: layrb(:,:)
      character (len=txtlen) :: string
      character (len=txtlen), pointer :: laytxt(:)
      intent (in) :: txtlen
      intent (out) :: layer
      intent (inout) :: string, layanz, laytxt, layrb
      end subroutine ftlay
end interface


interface
      subroutine funkt(x,p)
      use femtypes
      implicit none
      integer (I4B) p
      complex (DPC) x(:)
      intent (in) :: x, p
      end subroutine funkt
end interface


interface
      subroutine genele(ez,ezalt,ele,nb,geb,kdim,bz,                    &
     &  pkt1,pkt2,pkt3,prlist,prz)
      use femtypes
      implicit none
      integer (I4B) ez,ezalt,ele(:,:),nb(:,:),geb(:),kdim,bz
      integer (I4B) pkt1,pkt2,pkt3,prlist(:,:),prz
      intent (in) :: pkt1, pkt2, pkt3, ezalt, kdim, bz
      intent (out) :: nb, geb
      intent (inout) :: ez, ele, prlist, prz
      end subroutine genele
end interface


interface
      subroutine geng1(startz,startd,gbz,gzz,bzi,bzip,zki,xbk,ybk,      &
     &  area1,area2,negae,plp,plist)
      use femtypes
      implicit none
      integer (I4B) :: startz, startd, gbz, gzz, zki(:,:), negae
      integer (I4B), pointer :: bzi(:), bzip(:), plp(:), plist(:)
      real (DP) :: xbk(:), ybk(:)
      real (DP), pointer :: area1(:), area2(:)
      intent (in) :: startz, startd, gzz, zki, xbk, ybk
      intent (out) :: bzi, area1, area2, plist
      intent (inout) :: gbz, bzip, negae, plp
      end subroutine geng1
end interface


interface
      subroutine geng2(gbz,bzi,bzip,bzil,zki,xbk,ybk,area1,area2,     &
     &  negae,plp,plist)
      use femtypes
      implicit none
      integer (I4B) :: gbz
      integer (I4B) :: bzi(:), bzip(:), bzil(:), zki(:,:)
      integer (I4B) :: negae, plp(:), plist(:)
      real (DP) :: xbk(:), ybk(:), area1(:), area2(:)
      intent (in) :: gbz, zki, xbk, ybk, area1
      intent (out) :: bzi, bzip, bzil, area2, plist, plp
      intent (inout) :: negae
      end subroutine geng2
end interface


interface
      subroutine gengeb(gbz,gzz,gkz,bzi,bzip,bzil,zki,xbk,ybk,area1)
      use femtypes
      implicit none
      integer (I4B) :: gbz, gzz ,gkz, maxgbz, zki(:,:)
      integer (I4B), pointer :: bzi(:), bzip(:), bzil(:)
      real (DP) :: xbk(:), ybk(:)
      real (DP), pointer :: area1(:)
      intent (in) :: gzz, gkz, zki, xbk, ybk
      intent (out) :: gbz, bzi, bzip, bzil, area1
      end subroutine gengeb
end interface


interface
      subroutine genkz(zpz, zpp, branchpt, brnodes, bzi, bzil, bzip, meshsize)
      use femtypes
      implicit none
      integer (I4B) :: zpz(:), bzi(:), bzil(:), bzip(:)
      real (DP) :: meshsize(:)
      integer (I4B), pointer :: branchpt(:), brnodes(:)
      real (DP) :: zpp(:)
      intent (in) :: zpp, bzi, bzil, bzip, meshsize
      intent (inout) :: zpz
      end subroutine genkz
end interface


interface
      subroutine gennet(gbz,gzz,gkz,bzi,bzip,bzil,zki,zpz,zpp,zrb,      &
     &  xbk,ybk,kzi,xk,yk,kz,kdim,plist,prlist,ele,ez,geb,nb,           &
     &  xmin,xmax,ymin,ymax,plist2,plp,plp2)
      use femtypes
      implicit none
      integer (I4B) gkz, gzz, gbz, bzi(:), bzip(:), bzil(:), zki(:,:)
      integer (I4B) zpz(:), zrb(:), kzi(:), kz, kdim, plist(:)
      integer (I4B) prlist(:,:), ele(:,:), ez, geb(:), nb(:,:)
      integer (I4B) plist2(:), plp(:), plp2(:)
      real (DP) zpp(:)
      real (DP) xbk(:),ybk(:),xk(:),yk(:),xmin,xmax,ymin,ymax
      intent (in) :: gbz, gzz, gkz, bzi, bzip, bzil, zki, zpz, zpp, zrb
      intent (in) :: xbk, ybk, kdim
      intent (out) :: kzi, xk, yk, kz, ele, ez, geb, nb, xmin, xmax
      intent (out) :: ymin, ymax, prlist, plist, plist2, plp, plp2
      end subroutine gennet
end interface


interface
      subroutine gennode(i,j,nodelist,liste5,nregen,fehlkn)
      use femtypes
      implicit none
      integer (I4B) i, j, nodelist(:,:), liste5(:), fehlkn
      logical nregen 
      intent (in) :: i, j
      intent (out) :: nodelist, liste5, nregen
      intent (inout) :: fehlkn
      end subroutine gennode
end interface


interface
      subroutine genrcm(neqns, xadj, adjncy, perm)
      use femtypes
      implicit none
      integer (I4B) neqns, xadj(:), adjncy(:), perm(:)
      intent (in) :: neqns, adjncy
      intent (out) :: perm
      intent (inout) :: xadj
      end subroutine genrcm
end interface


interface
      pure subroutine get1Dintegpoints(intorder, xtab, weig, npkt, errcode)
      use femtypes
      implicit none
      integer (I4B) intorder, npkt, errcode
      real (DP), allocatable :: xtab(:), weig(:)
      intent (in) :: intorder
      intent (out) :: npkt, errcode, xtab, weig
      end subroutine get1Dintegpoints
end interface


interface
      pure subroutine get1Dintpolpoints(norder, xtab, errcode)
      use femtypes
      implicit none
      integer (I4B) norder, errcode
      real (DP) :: xtab(:)
      intent (in) :: norder
      intent (out) :: xtab, errcode
      end subroutine get1Dintpolpoints
end interface


interface
      pure subroutine get2Dintegpoints(intorder,npkt,weig,lambda,errcode)
      use femtypes
      implicit none
      integer (I4B) intorder, npkt, errcode
      real (DP), allocatable :: weig(:), lambda(:,:)
      intent (in) :: intorder
      intent (out) :: npkt, errcode, weig, lambda
      end subroutine get2Dintegpoints
end interface


interface
      subroutine getbcval2D( branch, nature, xs, ys, pval, qval, elem )
      use femtypes
      implicit none
      integer (I4B) :: branch, nature
      integer (I4B) , optional :: elem
      real (DP) :: xs, ys
      complex (DPC) :: pval
      complex (DPC), optional :: qval(:)
      intent (in)   :: branch, nature, xs, ys, elem
      intent (out)  :: pval, qval
      end subroutine getbcval2D
end interface


interface
      subroutine getfile(file,path)
      use femtypes
      implicit none
      character (len=*) file, path
      intent (out) file, path
      end subroutine getfile
end interface


interface
      subroutine getmat(matname,imat)
      use femtypes
      implicit none
      integer (I4B) imat
      character(len=*) matname 
      intent (in) :: matname
      intent (out) :: imat
      end subroutine getmat
end interface


interface
      subroutine getnue(i,x,xn,yn,e,nuex,nuey,dnyb2x,dnyb2y,            &
     &  bf,nuef,y2nues,n,ibegin)
      use femtypes
      implicit none
      integer (I4B) i,n,ibegin,e(:,:)
      real (SP) bf(:),nuef(:),y2nues(:)
      real (DP) nuex,nuey,dnyb2x,dnyb2y,xn(:),yn(:)
      complex (DPC)  x(:)
      intent (in) :: i,x,xn,yn,e,bf,nuef,y2nues,n,ibegin
      intent (out) :: nuex,nuey,dnyb2x,dnyb2y
      end subroutine getnue
end interface


interface getpostsetting
      subroutine getpostsetting_c(key,value)
      use femtypes
      implicit none
      character (len=*) key
      character (len=*) value
      intent (in) key
      intent (out) value
      end subroutine getpostsetting_c

      subroutine getpostsetting_i(key,value)
      use femtypes
      implicit none
      character (len=*) key
      integer (I4B) value
      intent (in) key
      intent (out) value
      end subroutine getpostsetting_i

      subroutine getpostsetting_f(key,value)
      use femtypes
      implicit none
      character (len=*) key
      real (DP) value
      intent (in) key
      intent (out) value
      end subroutine getpostsetting_f
end interface


interface getsetting
      recursive subroutine getsetting_c(key,value)
      use femtypes
      implicit none
      character (len=*) key
      character (len=*) value
      intent (in) key
      intent (out) value
      end subroutine getsetting_c

      recursive subroutine getsetting_i(key,value)
      use femtypes
      implicit none
      character (len=*) key
      integer (I4B) value
      intent (in) key
      intent (out) value
      end subroutine getsetting_i

      recursive subroutine getsetting_f(key,value)
      use femtypes
      implicit none
      character (len=*) key
      real (DP) value
      intent (in) key
      intent (out) value
      end subroutine getsetting_f
end interface


interface
      subroutine getsettingfile(settingfile)
      implicit none
      character (len=*) :: settingfile
      intent (out) :: settingfile
      end subroutine getsettingfile
end interface


interface
      subroutine getxyz(x,y,z)
      use femtypes
      implicit none
      real (DP) x,y,z
      intent (out) :: x, y, z
      end subroutine getxyz
end interface


interface
      subroutine glaett(n,e,xn,yn,en,geb,kzi,p,nregen,liste5,           &
     &  fehlkn,zrb,zki,xbk,ybk,loops)
      use femtypes
      implicit none
      integer (I4B) n, e(:,:), en(:,:), geb(:), kzi(:), p
      integer (I4B) liste5(:), fehlkn, zrb(:,:), zki(:,:), loops
      real (DP) xn(:), yn(:), xbk(:), ybk(:)
      logical nregen
      intent (in) :: n, geb, kzi, p, zrb, zki, xbk, ybk, loops
      intent (inout) :: e, xn, yn, en, nregen, liste5, fehlkn
      end subroutine glaett
end interface


interface
      pure subroutine gradshape(l, gl, polylo, polyhi, vec)
      use femtypes
      implicit none
      integer (I4B) polylo, polyhi
      real (DP) l(:), gl(:)
      real (DP) vec(:)
      intent (in) :: l, gl, polylo, polyhi
      intent (out) :: vec
      end subroutine gradshape
end interface


interface
      subroutine graf(a,nc,phi,amax,amin,h,x,y,e,ie,p,en)
      use femtypes
      implicit none
      complex (DPC) a(:)
      real (DP) h,x(:),y(:),phi,amax,amin
      integer (I4B) nc,e(:,:),ie,p,en(:,:)
      intent (in):: a,nc,phi,amax,amin,h,x,y,e,ie,p,en
      end subroutine graf
end interface


interface
      subroutine hadapt(ext,int,ende,error,ende2,epsgl)
      use femtypes
      implicit none
      real (DP) error, epsgl
      logical ende, ende2, ext,int
      intent (in) :: ext, int
      intent (out) :: ende, error, ende2, epsgl
      end subroutine hadapt
end interface


interface
      subroutine header(accur,xmin,xmax,ymin,ymax,unitfac)
      use femtypes
      implicit none
      real (DP) :: accur, xmax, xmin, ymax, ymin, unitfac
      intent (out) :: accur, xmax, xmin, ymax, ymin, unitfac
      end subroutine header
end interface


interface
      subroutine hfspkt(xp,yp,x1,y1,x2,y2,xmp,ymp,flaec2,length,ok)
      use femtypes
      implicit none
      real (DP) xp,yp,x1,y1,x2,y2,xmp,ymp,flaec2,length
      logical ok
      intent (in) :: xp, yp, x1, y1, x2, y2
      intent (out) :: xmp, ymp, flaec2, length, ok
      end subroutine hfspkt
end interface


interface
      pure subroutine hi2low(string,length)
      use femtypes
      implicit none
      integer (I4B) :: length
      character (len=*) :: string
      intent (in) :: length
      intent (inout) :: string
      end subroutine hi2low
end interface


interface
      subroutine hline(i,e,x,y,z,zcolor,ncolor,proj,transf)
      use femtypes
      implicit none
      real (DP) x(:), y(:), z(:) ,zcolor(:)
      real (DP), optional ::  transf(4,4)
      integer (I4B) i, e(:,:), ncolor
      logical proj
      intent (in):: i, e, x, y, z, zcolor, ncolor, transf, proj
      end subroutine hline
end interface


interface
      subroutine hpadapt(ext,int,ende,error,ende2,epsgl)
      use femtypes
      implicit none
      real (DP) error, epsgl
      logical ende, ende2, ext, int
      intent (in) :: ext, int
      intent (out) :: ende, error, ende2, epsgl
      end subroutine hpadapt
end interface


interface
      subroutine hpadapt_heuveline(ext,int,ende,error,ende2,epsgl)
      use femtypes
      implicit none
      real (DP) error, epsgl
      logical ende, ende2, ext, int
      intent (in) :: ext, int
      intent (out) :: ende, error, ende2, epsgl
      end subroutine hpadapt_heuveline
end interface


interface
      subroutine hpadapt_melenk(ext,int,ende,error,ende2,epsgl)
      use femtypes
      implicit none
      real (DP) error, epsgl
      logical ende, ende2, ext, int
      intent (in) :: ext, int
      intent (out) :: ende, error, ende2, epsgl
      end subroutine hpadapt_melenk
end interface


interface
      subroutine hpadapt_top5(ext,int,ende,error,ende2,epsgl)
      use femtypes
      implicit none
      real (DP) error, epsgl
      logical ende, ende2, ext, int
      intent (in) :: ext, int
      intent (out) :: ende, error, ende2, epsgl
      end subroutine hpadapt_top5
end interface


interface
      subroutine hpkt(j,ie,e,xn,yn,xmp,ymp,flaech2,weight)
      use femtypes
      implicit none
      real (DP) xn(:), yn(:), flaech2, xmp, ymp, weight
      integer (I4B) j, ie, e(:,:)
      intent (in) :: j, ie, e, xn, yn
      intent (out) ::  xmp, ymp, flaech2, weight
      end subroutine hpkt
end interface


interface
      subroutine hlinexy(i,e,x,y,z,z0,delta,proj,direction,transf)
      use femtypes
      implicit none
      real (DP) x(:), y(:), z(:), z0, delta
      real (DP), optional ::  transf(4,4)
      integer (I4B) i, e(:,:), direction
      logical proj
      intent (in):: i, e, x, y, z, z0, delta, transf, proj
      end subroutine hlinexy
end interface


interface
      subroutine hsb_to_rgb(hue,saturation,brightness,rgb)
      use femtypes
      implicit none
      real (DP) hue, saturation, brightness
      integer (I4B) rgb
      intent (in) :: hue, saturation, brightness
      intent (out) :: rgb
      end subroutine hsb_to_rgb
end interface


interface
      subroutine hunt(xx,x,jlo)
      use femtypes
      implicit none
      integer (I4B) :: jlo
      real (DP) :: x, xx(:)
      intent (in) :: x, xx
      intent (inout) :: jlo
      end subroutine hunt
end interface


interface
      subroutine idbvip(md,ndp,xd,yd,zd,nip,xi,yi,zi)
      implicit none
      integer :: md, ndp, nip
      real ( kind = 8 ) :: xd(ndp), yd(ndp), zd(ndp)
      real ( kind = 8 ) :: xi(nip), yi(nip), zi(nip)
      intent (in) :: md, ndp, xd, yd, zd, nip, xi, yi
      intent (out) :: zi
      end subroutine idbvip
end interface


interface
      subroutine indxf(anzknt,anzzwg,xpoint,ypoint,anzgeb,xgeb,ygeb,    &
     &  txtgeb,zki,zpz,lzrb,length,nohead,xmin,xmax,ymin,ymax,          &
     &  lytabl,accur,scalfk,lplayt,layanz,laytxt,layrb,snapof,txtlen)
      use femtypes
      implicit none
      integer (I4B) :: anzknt, anzzwg, anzgeb, layanz, snapof, txtlen
      integer (I4B), pointer :: zki(:,:), layrb(:,:), lplayt(:), zpz(:)
      integer (I4B), pointer :: lzrb(:)
      real (DP) :: xmin, xmax, ymin, ymax, accur, fakt1, scalfk
      real (DP), pointer :: xgeb(:), ygeb(:), xpoint(:), ypoint(:), length(:)
      logical :: nohead, lytabl
      character (len=txtlen), pointer :: txtgeb(:), laytxt(:)
      intent (in) :: txtlen
      intent (out) :: anzknt, anzzwg, xpoint, ypoint, anzgeb, zki
      intent (out) :: zpz, lzrb, length, nohead, xmin
      intent (out) :: xmax, ymin, ymax, lytabl, accur, scalfk, lplayt
      intent (out) :: layanz, laytxt, layrb, snapof
      intent (inout) :: txtgeb, xgeb, ygeb
      end subroutine indxf
end interface


interface
      subroutine initialize()
      use femtypes
      implicit none
      end subroutine initialize
end interface


interface
      subroutine inkrp(x1,y1,x2,y2,x3,y3,xpi,ypi,flaech2,rho,           &
     &                 xpa,ypa,rausen)
      use femtypes
      implicit none
      real (DP) x1,y1,x2,y2,x3,y3,xpi,ypi,flaech2,rho,xpa,ypa,rausen
      intent (in):: x1,y1,x2,y2,x3,y3
      intent (out):: xpi,ypi,flaech2,rho,xpa,ypa,rausen
      end subroutine inkrp
end interface


interface
      logical function innen(xk,yk,plist,pz,num1,num2,num3)
      use femtypes
      implicit none
      real (DP) xk(:),yk(:)
      integer (I4B) plist(:),pz
      integer (I4B) num1,num2,num3
      intent (in) :: xk,yk,plist,pz,num1,num2,num3
      end function innen
end interface


interface
      subroutine intbridges(head,nbds,numnodes)
      use femtypes
      use triangtypes
      implicit none
      integer (I4B)  nbds
      integer (I4B) :: numnodes(:)
      type (arrpt), pointer :: head(:)
      intent (in) :: nbds
      intent (inout) :: numnodes
      end subroutine intbridges
end interface


interface invertmat
      pure subroutine invertmat_dpc(matrix,inverted)
      use femtypes
      implicit none
      complex (DPC), intent (in) :: matrix(2,2)
      complex (DPC), intent (out) :: inverted(2,2)
      end subroutine invertmat_dpc

      pure subroutine invertmat_r(matrix,inverted)
      use femtypes
      implicit none
      real (DP), intent (in) :: matrix(2,2)
      real (DP), intent (out) :: inverted(2,2)
      end subroutine invertmat_r
end interface


interface invertmat3
      pure subroutine invertmat3_dpc(matrix,inverted)
      use femtypes
      implicit none
      complex (DPC), intent (in) :: matrix(3,3)
      complex (DPC), intent (out) :: inverted(3,3)
      end subroutine invertmat3_dpc

      pure subroutine invertmat3_r(matrix,inverted)
      use femtypes
      implicit none
      real (DP), intent (in) :: matrix(3,3)
      real (DP), intent (out) :: inverted(3,3)
      end subroutine invertmat3_r
end interface


interface
      subroutine ipop(ip,istack,lst,pttop,ptend,undrfl)
      use femtypes
      implicit none
      integer (I4B) ip,istack(:),lst,pttop,ptend
      logical undrfl
      intent (in) :: lst, ptend
      intent (out) :: ip, undrfl
      intent (inout) :: istack, pttop
      end subroutine ipop
end interface


interface
      subroutine ipush(ip,istack,lst,pttop,ptend,overfl)
      use femtypes
      implicit none
      integer (I4B) ip,istack(:),lst,pttop,ptend
      logical overfl
      intent (in) :: ip, lst, pttop
      intent (out) :: overfl
      intent (inout) :: istack, ptend
      end subroutine ipush
end interface


interface
      subroutine itrhstout(nnzero,iter,resvec,ndof)
      use femtypes
      implicit none
      real (DP) resvec(:)
      integer (I4B) :: nnzero, iter, ndof
      intent (in) :: nnzero, iter, resvec, ndof
      end subroutine itrhstout
end interface


interface
      subroutine juld2cal(jd,year,mon,day,month,dayofweek,hour,minute,second)
      use femtypes
      implicit none
      real (DP) jd
      integer (I4B) year, mon, day, hour, minute, second
      character(len=*) dayofweek, month
      intent (in) jd
      intent (out) year, mon, day, month, dayofweek, hour, minute, second
      end subroutine juld2cal
end interface


interface
      subroutine kern3(pliste,start,pzahl,xn,yn,xk,yk,pz)
      use femtypes
      implicit none
      integer (I4B) pliste(:), pzahl, pz ,start
      real (DP) xn(:), yn(:), xk(:), yk(:)
      intent (in) :: pliste, start, pzahl, xn, yn
      intent (out) :: xk, yk, pz
      end subroutine kern3
end interface


interface
      subroutine knsort(xn,yn,x,kzi,e,n,p,liste)
      use femtypes
      implicit none
      integer (I4B) liste(:),n,p,kzi(:),e(:,:)
      real (DP) xn(:),yn(:)
      complex (DPC) x(:)
      intent (in) :: n, p, liste
      intent (inout) :: xn, yn, x, kzi, e
      end subroutine knsort
end interface


interface
      subroutine kreipf(px,py,xp,yp,maxlen,vx,vy,acolor)
      use femtypes
      implicit none
      integer (I4B) px, py, acolor
      real (DP) maxlen
      real (DP) xp(:,:), yp(:,:), vx(:,:), vy(:,:)
      intent (in) :: px, py, xp, yp, maxlen, vx, vy, acolor
      end subroutine kreipf
end interface


interface
      pure subroutine lam2xy(lambda,elem,xt,yt,xn,yn,e)
      use femtypes
      implicit none
      integer (I4B) elem, e(:,:)
      real (DP) xt, yt, lambda(3), xn(:), yn(:)
      intent (in) :: lambda, xn, yn, elem, e
      intent (out) :: xt, yt
      end subroutine lam2xy
end interface


interface
      subroutine laytab(layanz,laytxt,layrb,txtlen)
      use femtypes
      implicit none
      integer (I4B) :: layanz, txtlen
      integer (I4B), pointer :: layrb(:,:)
      character (len=txtlen), pointer :: laytxt(:)
      intent (in) :: txtlen
      intent (out) :: layanz, laytxt, layrb
      end subroutine laytab
end interface


interface lcsr2coo
      subroutine lcsr2coo_DPC(diag, lower, upper, symm, lia, lja, ao, ir, jc)
      use femtypes
      implicit none
      complex (DPC), pointer:: diag(:), lower(:), upper(:), ao(:)
      integer (I4B), pointer :: lia(:), lja(:), ir(:), jc(:)
      logical symm
      intent(in) :: symm
      intent (out) :: ao, ir, jc
      intent (inout) :: diag, lower, upper, lia, lja
      end subroutine lcsr2coo_DPC

      subroutine lcsr2coo_SP(diag, lower, upper, symm, lia, lja, ao, ir, jc)
      use femtypes
      implicit none
      real (SP), pointer:: diag(:), lower(:), upper(:), ao(:)
      integer (I4B), pointer :: lia(:), lja(:), ir(:), jc(:)
      logical symm
      intent(in) :: symm
      intent (out) :: ao, ir, jc
      intent (inout) :: diag, lower, upper, lia, lja
      end subroutine lcsr2coo_SP

      subroutine lcsr2coo_DP(diag, lower, upper, symm, lia, lja, ao, ir, jc)
      use femtypes
      implicit none
      real (DP), pointer:: diag(:), lower(:), upper(:), ao(:)
      integer (I4B), pointer :: lia(:), lja(:), ir(:), jc(:)
      logical symm
      intent(in) :: symm
      intent (out) :: ao, ir, jc
      intent (inout) :: diag, lower, upper, lia, lja
      end subroutine lcsr2coo_DP
end interface


interface lcsr2csc
      subroutine lcsr2csc_DPC(diag, lower, upper, symm, lia, lja, acsc, ia, ja)
      use femtypes
      implicit none
      complex(DPC), pointer :: diag(:), lower(:), upper(:), acsc(:)
      integer(I4B), pointer :: lia(:), lja(:), ja(:), ia(:)
      logical symm
      intent(in) :: symm
      intent(out) :: acsc, ia, ja
      intent(inout) :: diag, lower, upper, lia, lja
      end subroutine lcsr2csc_DPC

      subroutine lcsr2csc_DP(diag, lower, upper, symm, lia, lja, acsc, ia, ja)
      use femtypes
      implicit none
      real(DP), pointer :: diag(:), lower(:), upper(:), acsc(:)
      integer(I4B), pointer :: lia(:), lja(:), ja(:), ia(:)
      logical symm
      intent(in) :: symm
      intent(out) :: acsc, ia, ja
      intent(inout) :: diag, lower, upper, lia, lja
      end subroutine lcsr2csc_DP
end interface


interface lcsr2csr
      subroutine lcsr2csr_DP(diag, lower, upper, symm, lia, lja, acsr, ia, ja)
      use femtypes
      implicit none
      real(DP), pointer :: diag(:), lower(:), upper(:), acsr(:)
      integer(I4B), pointer :: lia(:), lja(:), ja(:), ia(:)
      logical symm
      intent(in) :: symm
      intent(out) :: acsr, ia, ja
      intent(inout) :: diag, lower, upper, lia, lja
      end subroutine lcsr2csr_DP

      subroutine lcsr2csr_DPC(diag, lower, upper, symm, lia, lja, acsr, ia, ja)
      use femtypes
      implicit none
      complex(DPC), pointer :: diag(:), lower(:), upper(:), acsr(:)
      integer(I4B), pointer :: lia(:), lja(:), ja(:), ia(:)
      logical symm
      intent(in) :: symm
      intent(out) :: acsr, ia, ja
      intent(inout) :: diag, lower, upper, lia, lja
      end subroutine lcsr2csr_DPC
end interface


interface lcsr2lcsr
      subroutine lcsr2csr(diag, lower, upper, lia, lja, acsr, ia, ja)
      use femtypes
      implicit none
      complex(DPC), pointer :: diag(:), lower(:), upper(:), acsr(:)
      integer(I4B), pointer :: lia(:), lja(:), ja(:), ia(:)
      intent(out) :: acsr, ia, ja
      intent(inout) :: diag, lower, upper, lia, lja
      end subroutine lcsr2csr

      subroutine lcsr2lcsr_ZZ(diag_in, lower_in, upper_in, rhs_in, symm, diag_out, lower_out, upper_out, rhs_out)
      use femtypes
      implicit none
      complex(DPC), pointer :: diag_in(:), lower_in(:), upper_in(:), rhs_in(:)
      complex(DPC), pointer :: diag_out(:), lower_out(:), upper_out(:), rhs_out(:)
      logical symm
      intent(in) :: symm
      intent(out) :: diag_out, lower_out, upper_out, rhs_out
      intent(inout) :: diag_in, lower_in, upper_in, rhs_in
      end subroutine lcsr2lcsr_ZZ

      subroutine lcsr2lcsr_ZD(diag_in, lower_in, upper_in, rhs_in, symm, diag_out, lower_out, upper_out, rhs_out)
      use femtypes
      implicit none
      complex(DPC), pointer :: diag_in(:), lower_in(:), upper_in(:), rhs_in(:)
      real(DP), pointer :: diag_out(:), lower_out(:), upper_out(:), rhs_out(:)
      logical symm
      intent(in) :: symm
      intent(out) :: diag_out, lower_out, upper_out, rhs_out
      intent(inout) :: diag_in, lower_in, upper_in, rhs_in
      end subroutine lcsr2lcsr_ZD

      subroutine lcsr2lcsr_ZC(diag_in, lower_in, upper_in, rhs_in, symm, diag_out, lower_out, upper_out, rhs_out)
      use femtypes
      implicit none
      complex(DPC), pointer :: diag_in(:), lower_in(:), upper_in(:), rhs_in(:)
      complex(SPC), pointer :: diag_out(:), lower_out(:), upper_out(:), rhs_out(:)
      logical symm
      intent(in) :: symm
      intent(out) :: diag_out, lower_out, upper_out, rhs_out
      intent(inout) :: diag_in, lower_in, upper_in, rhs_in
      end subroutine lcsr2lcsr_ZC

      subroutine lcsr2lcsr_ZS(diag_in, lower_in, upper_in, rhs_in, symm, diag_out, lower_out, upper_out, rhs_out)
      use femtypes
      implicit none
      complex(DPC), pointer :: diag_in(:), lower_in(:), upper_in(:), rhs_in(:)
      real(SPC), pointer :: diag_out(:), lower_out(:), upper_out(:), rhs_out(:)
      logical symm
      intent(in) :: symm
      intent(out) :: diag_out, lower_out, upper_out, rhs_out
      intent(inout) :: diag_in, lower_in, upper_in, rhs_in
      end subroutine lcsr2lcsr_ZS
end interface


interface
      subroutine lcsrmatrhsout(n,ia,ja,diag,lower,upper,rhs1,ascii)
      use femtypes
      implicit none
      integer (I4B) n
      complex (DPC) rhs1(:), diag(:) ,lower(:) ,upper(:)
      integer (I4B) ia(:), ja(:)
      logical ascii
      intent (in) n, ia, ja, diag, lower, upper, rhs1, ascii
      end subroutine lcsrmatrhsout
end interface


interface lcsrzeroremover
      subroutine lcsrzeroremover_ZZ(diag, lower, upper, lia, lja, symm, method, eps, rhs, x, eps_x)
      use femtypes
      implicit none
      complex(DPC), pointer :: diag(:), lower(:), upper(:)
      complex(DPC), pointer, optional :: x(:), rhs(:)
      real(DP), optional :: eps, eps_x
      integer(I4B), pointer :: lia(:), lja(:)
      integer (I4B) :: method
      logical :: symm
      intent(in) :: symm, method, eps_x, x
      intent(inout) :: diag, lower, upper, lia, lja, rhs
      end subroutine lcsrzeroremover_ZZ

      subroutine lcsrzeroremover_CC(diag, lower, upper, lia, lja, symm, method, eps, rhs, x, eps_x)
      use femtypes
      implicit none
      complex(SPC), pointer :: diag(:), lower(:), upper(:)
      complex(SPC), pointer, optional :: x(:), rhs(:)
      real(DP), optional :: eps, eps_x
      integer(I4B), pointer :: lia(:), lja(:)
      integer (I4B) :: method
      logical :: symm
      intent(in) :: symm, method, eps_x, x
      intent(inout) :: diag, lower, upper, lia, lja, rhs
      end subroutine lcsrzeroremover_CC

      subroutine lcsrzeroremover_DD(diag, lower, upper, lia, lja, symm, method, eps, rhs, x, eps_x)
      use femtypes
      implicit none
      real(DP), pointer :: diag(:), lower(:), upper(:)
      real(DP), pointer, optional :: x(:), rhs(:)
      real(DP), optional :: eps, eps_x
      integer(I4B), pointer :: lia(:), lja(:)
      integer (I4B) :: method
      logical :: symm
      intent(in) :: symm, method, eps_x, x
      intent(inout) :: diag, lower, upper, lia, lja, rhs
      end subroutine lcsrzeroremover_DD

      subroutine lcsrzeroremover_SS(diag, lower, upper, lia, lja, symm, method, eps, rhs, x, eps_x)
      use femtypes
      implicit none
      real(SP), pointer :: diag(:), lower(:), upper(:)
      real(SP), pointer, optional :: x(:), rhs(:)
      real(DP), optional :: eps, eps_x
      integer(I4B), pointer :: lia(:), lja(:)
      integer (I4B) :: method
      logical :: symm
      intent(in) :: symm, method
      intent(inout) :: diag, lower, upper, lia, lja
      end subroutine lcsrzeroremover_SS
end interface


interface
      subroutine leftpolygon(xa,ya,xb,yb,pin,xi,yi,pout,xo,yo)
      use femtypes
      implicit none
      integer (I4B) pin, pout
      real (DP) xa, ya, xb, yb, xi(:), yi(:), xo(:), yo(:)
      intent (in) :: xa, ya, xb, yb, pin, xi, yi
      intent (out) :: pout, xo, yo
      end subroutine leftpolygon
end interface


interface
      subroutine lese(bzi,bzip,bzil,zpz,zpp,matname,krb,meshsize,       &
     &                lzrb,layanz,layrb,laytxt,lalrb,lbtrb)
      use femtypes
      implicit none
      integer (I4B) :: layanz
      integer (I4B), pointer :: bzi(:), bzip(:), bzil(:), krb(:,:)
      integer (I4B), pointer :: zpz(:), lzrb(:), layrb(:,:)
      real (DP), pointer :: zpp(:), meshsize(:)
      complex (DPC), pointer :: lalrb(:,:), lbtrb(:,:)
      character (len=*), pointer :: matname(:), laytxt(:)
      intent (out) :: bzi, bzip, bzil, zpz, zpp, matname, krb, meshsize
      intent (out) :: lzrb, layanz, layrb, laytxt, lalrb, lbtrb
      end subroutine lese
end interface


interface
      subroutine lin(ok, jdmesh, jdsolu, gilt, eleinfo, n, ndof, omega, &
     &               epsgl, nnat, ep, eg, x, fem_accuracy)
      use femtypes
      implicit none
      integer (I4B) :: ndof, n, nnat
      integer (I4B), pointer :: ep(:,:)
      complex (DPC), pointer :: x(:)
      type (ARRPTRI), pointer :: eg(:,:)
      real (DP) :: omega, epsgl, jdmesh, jdsolu, fem_accuracy
      logical ok, gilt, eleinfo
      intent (in) :: jdmesh, n, nnat
      intent (out) :: ok, jdsolu, gilt, eleinfo, ndof, omega
      intent (out) :: epsgl, fem_accuracy
      end subroutine lin
end interface


interface
      subroutine line(accur,fakt1,scalfk,laytxt,layrb,layanz,xpoint,    &
     &  ypoint,zki,lzrb,zpz,length,anzzwg,anzknt,ok,txtlen)
      use femtypes
      implicit none
      integer (I4B) :: anzzwg, anzknt, txtlen, layanz
      integer (I4B), pointer :: lzrb(:), zki(:,:), zpz(:), layrb(:,:)
      real (DP) :: accur, fakt1, scalfk
      real (DP), pointer :: xpoint(:), ypoint(:), length(:)
      logical :: ok
      character (len=txtlen), pointer :: laytxt(:)
      intent (in) :: accur, fakt1, scalfk, txtlen
      intent (out) :: xpoint, ypoint, zki, lzrb, zpz, length, ok
      intent (inout) :: laytxt, layrb, layanz, anzzwg, anzknt
      end subroutine line
end interface


interface
      subroutine linear_sparse_solve(numdof, lower_DPC, upper_DPC, diag_DPC, rhs_DPC, x, ia, ja, epsgl, resgl)
      use femtypes
      implicit none
      integer(I4B) :: numdof
      integer(I4B), pointer :: ia(:), ja(:)
      real(DP) :: epsgl, resgl
      complex(DPC), pointer :: lower_DPC(:), upper_DPC(:), diag_DPC(:), rhs_DPC(:), x(:)
      intent(in) :: numdof
      intent(out) :: epsgl, resgl
      intent(inout) :: lower_DPC, upper_DPC, diag_DPC, rhs_DPC, ia, ja, x
      end subroutine linear_sparse_solve
end interface


interface
      subroutine linfluid(ok, jdmesh, jdsolu, gilt, eleinfo, n, ndof, omega, &
     &  epsgl, ep, eg, xfluid, fem_accuracy)
      use femtypes
      implicit none
      integer (I4B) ndof, n
      integer (I4B), pointer :: ep(:,:)
      complex (DPC), pointer :: xfluid(:,:)
      type (ARRPTRI), pointer :: eg(:,:)
      real (DP) omega, epsgl, jdmesh, jdsolu, fem_accuracy
      logical ok, gilt, eleinfo
      intent (in) :: jdmesh, n
      intent (out) :: ok, jdsolu, gilt, eleinfo, ndof, omega
      intent (out) :: epsgl, fem_accuracy
      end subroutine linfluid
end interface


interface
      subroutine linstokes(ok, jdmesh, jdsolu, gilt, eleinfo, n, ndof, omega, &
     &  epsgl, ep, eg, xfluid, fem_accuracy)
      use femtypes
      implicit none
      integer (I4B) ndof, n !, npdof
      integer (I4B), pointer :: ep(:,:)
      complex (DPC), pointer :: xfluid(:,:)
      type (ARRPTRI), pointer :: eg(:,:)
      real (DP) omega, epsgl, jdmesh, jdsolu, fem_accuracy
      logical ok, gilt, eleinfo
      intent (in) :: jdmesh, n
      intent (out) :: ok, jdsolu, gilt, eleinfo, ndof, omega !, npdof
      intent (out) :: epsgl, fem_accuracy
      end subroutine linstokes
end interface


interface locate
      function locate_DP(xx,x)
      use femtypes
      implicit none
      real (DP) xx(:),x
      integer (I4B) locate_DP
      intent (in) :: x,xx
      end function locate_DP

      function locate_idx_DP(xx,x,indx)
      use femtypes
      implicit none
      real (DP) xx(:),x
      integer (I4B) locate_idx_DP, indx(:)
      intent (in) :: x,xx,indx
      end function locate_idx_DP
end interface


interface
      subroutine loesch(negae,merk,plp,plist,area2)
      use femtypes
      implicit none
      integer (I4B) :: negae, merk, plp(:), plist(:)
      real (DP) :: area2(:)
      intent (in) :: merk
      intent (inout) :: negae, plist, plp, area2
      end subroutine loesch
end interface


interface
      subroutine lokmat(gebiet,matzif,materi)
      use femtypes
      implicit none
      integer (I4B) gebiet,matzif(:),materi
      intent (in) :: gebiet,matzif
      intent (out) :: materi
      end subroutine lokmat
end interface


interface
      subroutine lout(gilt, eleinfo, n, ndof, omega, epsgl, nnat, ep,     &
     &           eg, x, fem_accuracy)
      use femtypes
      implicit none
      integer (I4B) :: n
      integer (I4B) :: ndof, nnat
      integer (I4B), optional :: ep(:,:)
      complex (DPC), optional :: x(:)
      type (ARRPTRI), optional :: eg(:,:)
      real (DP) :: omega, epsgl, fem_accuracy
      logical gilt, eleinfo
      intent (in) :: gilt, eleinfo, n, ndof, nnat, omega, epsgl, fem_accuracy
      end subroutine lout
end interface


interface
      pure subroutine low2hi(string,length)
      use femtypes
      implicit none
      integer (I4B) :: length
      character (len=*) :: string
      intent (in) :: length
      intent (inout) :: string
      end subroutine low2hi
end interface


interface
      subroutine lsmoth(kzi,p,xn,yn,e,n,en,itermin,nregen,liste5,       &
     &  fehlkn,zrb,zki,xbk,ybk)
      use femtypes
      implicit none
      integer (I4B) e(:,:), kzi(:), n, p, en(:,:), itermin
      integer (I4B) liste5(:), fehlkn, zrb(:,:), zki(:,:)
      real (DP) xn(:), yn(:), xbk(:), ybk(:)
      logical nregen
      intent (in) :: kzi, p, e, n, en, itermin, zrb, zki, xbk, ybk
      intent (out) :: nregen
      intent (inout) :: xn, yn, liste5, fehlkn
      end subroutine lsmoth
end interface


interface
      subroutine lswap(n,e,xn,yn,en,geb,talk,all,szahl,ic,slist)
      use femtypes
      implicit none
      integer (I4B) n, e(:,:), en(:,:), geb(:), talk, szahl, ic
      integer (I4B), optional, target :: slist(:)
      real (DP) xn(:), yn(:)
      logical all
      intent (in) :: n, xn, yn, geb, ic, talk, all
      intent (out) :: szahl
      intent (inout) :: e, en, slist
      end subroutine lswap
end interface


interface
      subroutine lswap1(e,en,xn,yn,elnum,ennum,z1,z2,swappd,ic)
      use femtypes
      implicit none
      integer (I4B) e(:,:), en(:,:), elnum, ennum, z1, z2, ic
      real (DP) xn(:), yn(:)
      logical swappd
      intent (in) :: xn, yn, elnum, ennum, z1, z2, ic
      intent (out) :: swappd
      intent (inout) :: e, en
      end subroutine lswap1
end interface


interface lubksb
      pure subroutine lubksb_c(a,indx,b)
      use femtypes
      implicit none
      complex (DPC) a(:,:), b(:)
      integer (I4B) indx(:)
      intent (in) :: a, indx
      intent (inout) :: b
      end subroutine lubksb_c

      pure subroutine lubksb_d(a,indx,b)
      use femtypes      
      implicit none      
      real (DP) :: a(:,:), b(:)
      integer (I4B) :: indx(:)
      intent (in) :: a, indx
      intent (inout) :: b
      end subroutine lubksb_d
end interface lubksb


interface ludcmp
      subroutine ludcmp_c(a,indx)
      use femtypes
      implicit none
      complex (DPC) a(:,:)
      integer (I4B) indx(:)
      intent (out) :: indx
      intent (inout) :: a
      end subroutine ludcmp_c

      subroutine ludcmp_d(a,indx)
      use femtypes
      implicit none
      real (DP) :: a(:,:)
      integer (I4B) indx(:)
      intent (out) :: indx
      intent (inout) :: a
      end subroutine ludcmp_d
end interface ludcmp


interface lusolver
      subroutine lusolver_c(a,b)
      use femtypes
      implicit none
      complex (DPC) :: a(:,:), b(:)
      intent (inout) :: a, b
      end subroutine lusolver_c

      subroutine lusolver_d(a,b)
      use femtypes
      implicit none
      real (DP) :: a(:,:), b(:) 
      intent (inout) :: a, b
      end subroutine lusolver_d
end interface lusolver


interface
      subroutine lwpolyline(accur,fakt1,scalfk,laytxt,layrb,layanz,     &
     &  xpoint,ypoint,zki,lzrb,zpz,length,anzzwg,anzknt,ok,txtlen)
      use femtypes
      implicit none
      integer (I4B) :: anzzwg, anzknt, txtlen, layanz
      integer (I4B), pointer :: lzrb(:), zki(:,:), zpz(:), layrb(:,:)
      real (DP) :: accur, fakt1, scalfk
      real (DP), pointer :: xpoint(:), ypoint(:), length(:)
      logical :: ok
      character (len=txtlen), pointer :: laytxt(:)
      intent (in) :: accur, fakt1, scalfk
      intent (out) :: xpoint, ypoint, zki, lzrb, zpz, length, ok
      intent (inout) :: laytxt, layrb, layanz, anzzwg, anzknt
      end subroutine lwpolyline
end interface


interface
      subroutine madjncy(e,en,n,xadj,adjncy)
      use femtypes
      implicit none
      integer (I4B)  e(:,:),en(:,:),xadj(:),n, adjncy(:)
      intent (in) :: e, en, n, xadj
      intent (out) :: adjncy
      end subroutine madjncy
end interface


interface mat2ten
      pure subroutine mat2ten_r4_DP(matrix, tensor, factors)
      use femtypes
      implicit none
      real (DP) :: matrix(6,6), tensor(3,3,3,3), factors(3)
      intent (in) :: matrix, factors
      intent (out) :: tensor
      end

      pure subroutine mat2ten_r4_DPC(matrix, tensor, factors)
      use femtypes
      implicit none
      real(DP) :: factors(3)
      complex(DPC) :: matrix(6,6), tensor(3,3,3,3)
      intent (in) :: matrix, factors
      intent (out) :: tensor
      end
end interface

 !TODO: Replace old tensor functions by new ones
interface mat2ten_r2
      pure function mat2ten_r2_dpc(mat)
      use femtypes
      implicit none
      complex (DPC) :: mat(6), mat2ten_r2_dpc(3,3)
      intent (in) :: mat
      end function mat2ten_r2_dpc

      pure function mat2ten_r2_dp(mat)
      use femtypes
      implicit none
      real (DP) :: mat(6), mat2ten_r2_dp(3,3)
      intent (in) :: mat
      end function mat2ten_r2_dp
end interface


interface mat2ten_r3
      pure function mat2ten_r3_dpc(mat)
      use femtypes
      implicit none
      complex (DPC) :: mat(6,3), mat2ten_r3_dpc(3,3,3)
      intent (in) :: mat
      end function mat2ten_r3_dpc

      pure function mat2ten_r3_dp(mat)
      use femtypes
      implicit none
      real (DP) :: mat(6,3), mat2ten_r3_dp(3,3,3)
      intent (in) :: mat
      end function mat2ten_r3_dp
end interface


interface mat2ten_r4
      pure function mat2ten_r4_dpc_old(mat)
      use femtypes
      implicit none
      complex (DPC) :: mat(6,6), mat2ten_r4_dpc_old(3,3,3,3)
      intent (in) :: mat
      end function mat2ten_r4_dpc_old

      pure function mat2ten_r4_dp_old(mat)
      use femtypes
      implicit none
      real (DP) :: mat(6,6), mat2ten_r4_dp_old(3,3,3,3)
      intent (in) :: mat
      end function mat2ten_r4_dp_old
end interface


interface
      subroutine matgeb(gbz,anzgeb,bzi,bzip,bzil,zki,area1,xpoint,      &
     &                  ypoint,xgeb,ygeb,txtgeb,matname,donetx)
      use femtypes
      implicit none
      integer (I4B) :: gbz, anzgeb, bzi(:), bzip(:), bzil(:)
      integer (I4B) :: zki(:,:), donetx(:)
      real (DP) :: xpoint(:), ypoint(:), xgeb(:), ygeb(:), area1(:)
      character (len=*) :: txtgeb(:), matname(:)
      intent (in) :: gbz, anzgeb, bzi, bzip, bzil, zki, area1
      intent (in) ::  xpoint, ypoint
      intent (out) :: matname
      intent (inout) :: txtgeb
      end subroutine matgeb
end interface


interface
      subroutine mathbout(n,ia,ja,csr)
      use femtypes
      implicit none
      integer (I4B) n
      complex (DPC) csr(:)
      integer (I4B) ia(:), ja(:)
      intent (in) n, ia, ja, csr
      end subroutine mathbout
end interface


interface
      subroutine matout(lower,upper,diag,b,ia,ja,p,n,                 &
     &  rlower,rupper,rdiag,rb,reell)
      use femtypes
      implicit none
      integer (I4B) ja(:),ia(:),p,n
      complex (SPC) lower(:),upper(:),diag(:),b(:)
      real (DP) rlower(:),rupper(:),rdiag(:),rb(:)
      logical reell
      intent (in) :: lower, upper, diag, b, ia, ja, p, n
      intent (in) :: rlower, rupper, rdiag, rb, reell
      end subroutine matout
end interface


interface
      subroutine matrhsout(n,ia,ja,csr,rhs1,ascii)
      use femtypes
      implicit none
      integer (I4B) n
      complex (DPC) rhs1(:), csr(:)
      integer (I4B) ia(:), ja(:)
      logical ascii
      intent (in) n, ia, ja, csr, rhs1, ascii
      end subroutine matrhsout
end interface


interface
      subroutine meshinput(ok,jd,reell,gilt,istlin,konfom,zyl,epsgl,resgl,resnl,error)
      use femtypes
      implicit none
      real (DP) epsgl, resgl, resnl, error, jd
      integer (I4B) ndof
      logical ok, reell, gilt, istlin, konfom, zyl

      intent (out) :: ok, reell, gilt, istlin, konfom, zyl, jd
      intent (out) :: epsgl, resgl, resnl, error
      end subroutine meshinput
end interface


interface
      subroutine minmesh(bzi,bzip,bzil,zpz,zpp, meshsize)
      use femtypes
      use triangtypes
      implicit none
      integer (I4B) bzi(:), bzip(:), bzil(:), zpz(:)
      real (DP) zpp(:), meshsize(:)
      intent (in) :: bzi, bzip, bzil, zpp, meshsize
      intent (inout) :: zpz
      end subroutine minmesh
end interface


!interface mminfo
!      subroutine  subroutine mminfo(iunit,rep,field,symm,rows,cols,nnz)
!      integer i, rows, cols, nnz, iunit
!      integer count
!      character mmhead*14
!      character mmtype*6
!      character rep*10
!      character field*7
!      character symm*19
!      character tmp1*1024
!      character tmp2*2
!end interface


!interface mmread
!      subroutine mmread(iunit,rep,field,symm,rows,cols,nnz,nnzmax, &
!                 &      indx,jndx,ival,rval,cval)
!      integer ival(*)
!      double precision rval(*)
!      complex cval(*)
!      double precision rpart,ipart
!      integer indx(*)
!      integer jndx(*)
!      integer i, rows, cols, nnz, nnzreq, nnzmax, iunit
!      integer count
!      character mmhead*15
!      character mmtype*6
!      character rep*10
!      character field*7
!      character symm*19
!      character tmp1*1024
!      character tmp2*2
!end interface


interface   
      subroutine mmwrite(ounit,rep,field,symm,rows,cols,nnz, &
                 &       indx,jndx,ival,rval,cval)
      use femtypes
      integer ival(*)
      double precision rval(*)
      complex (DPC) cval(*)
      integer indx(*)
      integer jndx(*)
      integer rows, cols, nnz, ounit
      character*(*)rep,field,symm
      end subroutine
end interface


interface mp
      subroutine mp_c(lower,upper,diag,b,x,ia,ja,n)
      use femtypes
      implicit none
      complex (SPC) lower(:),upper(:),diag(:),x(:),b(:)
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
      end subroutine mp_c

      subroutine mp_z(lower,upper,diag,b,x,ia,ja,n)
      use femtypes
      implicit none
      complex (DPC) lower(:),upper(:),diag(:),x(:),b(:)
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
      end subroutine mp_z

      subroutine mp_d(lower,upper,diag,b,x,ia,ja,n)
      use femtypes
      implicit none
      real (DP) lower(:),upper(:),diag(:),x(:),b(:)
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
      end subroutine mp_d

      subroutine mp_s(lower,upper,diag,b,x,ia,ja,n)
      use femtypes
      implicit none
      real (SP) lower(:),upper(:),diag(:),x(:),b(:)
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
      end subroutine mp_s
end interface


interface
      subroutine mplist(gkz,bz,bzi,bzip,bzil,zki,zpz,                   &
     &  plist,pz,plp,zahl)
      use femtypes
      implicit none
      integer (I4B) gkz,bz,bzi(:),bzip(:),bzil(:),zki(:,:)
      integer (I4B) zpz(:),plist(:),pz,plp(:),zahl
      intent (in) :: gkz, bz, bzi, bzip, bzil, zki, zpz
      intent (out) :: plist, plp, pz, zahl
      end subroutine mplist
end interface


interface
      subroutine mstruk(gbz,bzi,bzip,bzil,zki,zrb,xbk,ybk,              &
     &  xn,yn,p,e,n,xmin,xmax,ymin,ymax,krb,kzi)
      use femtypes
      implicit none
      integer (I4B) gbz, bzi(:), bzip(:), bzil(:), zki(:,:), zrb(:)
      integer (I4B) p, e(:,:), n, krb(:), kzi(:)
      real (DP) xn(:), yn(:), xbk(:), ybk(:)
      real (DP) xmin, xmax, ymin, ymax
      intent (in) :: gbz, bzi, bzip, bzil, zki, zrb, xbk, ybk
      intent (in) :: xn, yn, p, e, n, xmin, xmax, ymin, ymax, krb, kzi
      end subroutine mstruk
end interface


interface
      subroutine mtext(accur,scalfk,anzgeb,xgeb,ygeb,txtgeb,layanz,     &
     &  lplayt,laytxt,layrb,txtlen)
      use femtypes
      implicit none
      integer (I4B) :: anzgeb, txtlen, layanz
      integer (I4B), pointer :: layrb(:,:), lplayt(:)
      real (DP) :: accur, scalfk
      real (DP), pointer :: xgeb(:), ygeb(:)
      character (len=txtlen), pointer :: txtgeb(:), laytxt(:)
      intent (in) :: accur, scalfk, txtlen
      intent (out) :: lplayt
      intent (inout) :: xgeb, ygeb, txtgeb
      intent (inout) :: anzgeb, layanz, laytxt, layrb
      end subroutine mtext
end interface


interface multtensor
      pure subroutine multtensor_ECD(tensor1, tensor2, ten)
      use femtypes
      implicit none
      real (DP) :: tensor1(2,2,2,2), tensor2(2,2,2), ten(2,2,2)
      intent (in) :: tensor1, tensor2
      intent (out) :: ten
      end subroutine multtensor_ECD

      pure subroutine multtensor_EDC(tensor1, tensor2, ten)
      use femtypes
      implicit none
      real (DP) :: tensor1(2,2,2), tensor2(2,2,2,2), ten(2,2,2)
      intent (in) :: tensor1, tensor2
      intent (out) :: ten
      end subroutine multtensor_EDC

      pure subroutine multtensor_MCA(tensor1, tensor2, ten)
      use femtypes
      implicit none
      real (DP) :: tensor1(2,2,2,2), tensor2(2,2), ten(2,2)
      intent (in) :: tensor1, tensor2
      intent (out) :: ten
      end subroutine multtensor_MCA
end interface


interface
      subroutine mxadj(e,en,n,p,xadj)
      use femtypes
      implicit none
      integer (I4B)  e(:,:),en(:,:),xadj(:),n,p
      intent (in) :: e, en, n, p
      intent (out) :: xadj
      end subroutine mxadj
end interface


interface
      subroutine nachb(e,en,p,n)
      use femtypes
      implicit none
      integer (I4B) p, e(:,:), en(:,:), n
      intent (in) e, p, n
      intent (out) en
      end subroutine nachb
end interface


interface
      subroutine nachkn(p,n,kanz,kknot,e,en,kzi,zrb,zki,liste5)
      use femtypes
      implicit none
      integer (I4B) p, n, kanz(:), kknot(:), e(:,:), en(:,:), kzi(:)
      integer (I4B) zrb(:,:), zki(:,:), liste5(:)
      intent (in) :: p, n, e, en, kzi, zrb, zki, liste5
      intent (out) :: kanz, kknot
      end subroutine nachkn
end interface


interface
      subroutine netbew(xn,yn,e,rk1,n)
      use femtypes
      implicit none
      real (DP) xn(:),yn(:), rk1
      integer (I4B) e(:,:), n
      intent (in) :: xn, yn, e,  n
      intent (out) :: rk1
      end subroutine netbew
end interface


interface
      subroutine neukn(n1,n2,n3,n4,neuk,xn,yn,x,kzi,zki,p,xbk,ybk,      &
     &  geb1,geb2,arc,err)
      use femtypes
      implicit none
      integer (I4B) n1, n2, n3, n4, neuk, kzi(:), zki(:,:)
      integer (I4B) p,geb1,geb2, err
      real (DP) xn(:),yn(:),xbk(:),ybk(:)
      complex (DPC) x(:)
      logical arc
      intent (in) :: n1, n2, n3, n4, zki, xbk, ybk, geb1, geb2
      intent (out) :: neuk, arc, err
      intent (inout) :: xn, yn, x, kzi, p
      end subroutine neukn
end interface


interface
      subroutine newelements(ellist,nregen,fehlkn,liste5)
      use femtypes
      implicit none
      integer (I4B) :: fehlkn, ellist(:), liste5(:)
      logical nregen
      intent (in) ::  ellist
      intent (out) :: nregen, fehlkn, liste5
      end  subroutine newelements
end interface


interface
      subroutine newelementshp(father,ellist,nregen,fehlkn,liste5)
      use femtypes
      implicit none
      integer (I4B) :: fehlkn, ellist(:), liste5(:), father(:)
      logical nregen
      intent (in) ::  ellist
      intent (out) :: nregen, fehlkn, liste5
      intent (inout) :: father
      end  subroutine newelementshp
end interface


interface
      subroutine nfarbe(z,x,y,e,ie,palette,ncolor,zcolmin,zcolmax)
      use femtypes
      implicit none
      real (DP) x(:), y(:), z(:), zcolmin, zcolmax
      integer (I4B) e(:,:), ie, ncolor, palette(:)
      intent (in):: z, x, y, e, ie, palette, ncolor, zcolmin, zcolmax
      end subroutine nfarbe
end interface


interface
      subroutine nfarbep(palette,colorfill,ncolor,zcolmin,zcolmax,      &
     &  hlin,nlin,zlinmin,zlinmax,r,g,b,phi,fieldtype)
      use femtypes
      implicit none
      real (SP) :: r, g, b
      real (DP) :: zcolmin, zcolmax, zlinmin, zlinmax, phi
      integer (I4B) :: ncolor, nlin, palette(:)
      logical :: colorfill, hlin
      character (len=*) :: fieldtype
      intent (in):: palette, colorfill, ncolor, zcolmin, zcolmax, hlin, nlin
      intent (in):: zlinmin, zlinmax, r, g, b, phi, fieldtype
      end subroutine nfarbep
end interface


interface
      logical function ninnen(x1,y1,x2,y2,x3,y3,xk,yk)
      use femtypes
      implicit none
      real (DP) x1, y1, x2, y2, x3, y3, xk, yk
      intent (in) :: x1, y1, x2, y2, x3, y3, xk, yk
      end function ninnen
end interface


interface nonlinearsolver
      subroutine nonlinearsolver_newmark(acsr,x_at_t,csrvec,ia,ja,newmark_params,epsgl,resgl)
      use femtypes
      implicit none
      complex (DPC), pointer :: acsr(:), x_at_t(:), csrvec(:)
      integer (I4B), pointer :: ia(:), ja(:)
      real (DP) :: newmark_params(8),epsgl, resgl
      intent (in) :: newmark_params
      intent (out) :: epsgl, resgl
      end subroutine nonlinearsolver_newmark

      subroutine nonlinearsolver_stat(acsr,rhs,lower,upper,diag,ia,ja,epsgl,resgl)
      use femtypes
      implicit none
      complex (DPC), pointer :: acsr(:), rhs(:), lower(:), upper(:), diag(:)
      integer (I4B), pointer :: ia(:), ja(:)
      real (DP) :: epsgl, resgl
      intent (out) :: epsgl, resgl
      end subroutine nonlinearsolver_stat
end interface


interface
      subroutine nsurf(z,x,y,p,e,ie,palette,ncolor,zcolmin,zcolmax,     &
     &     zscale,color,hlin,cross,delcross,theta,phi)
      use femtypes
      implicit none
      real (DP) x(:), y(:), z(:), zcolmin, zcolmax, zscale, delcross 
      real (DP) phi, theta
      integer (I4B) p, e(:,:), ie, ncolor, palette(:)
      logical color, hlin, cross
      intent (in):: z, x, y, p, e, ie, palette, ncolor, zcolmin,         &
     &   zcolmax, zscale, color, hlin, cross, delcross, theta, phi
      end subroutine nsurf
end interface


interface
      subroutine numberpoly(dch,xdmin,xdmax,ydmin,ydmax,nature)
      use femtypes
      implicit none
      integer (I4B) :: nature
      real (DP) :: dch, xdmin, xdmax, ydmin, ydmax
      intent (in) :: dch, xdmin, xdmax, ydmin, ydmax, nature
      end subroutine numberpoly
end interface


interface
      subroutine objects
      use femtypes
      implicit none
      end subroutine objects
end interface


interface
      subroutine offcentre(x,y,alpha,minloc,betamax,xpa,ypa)
      use femtypes
      implicit none
      real (DP) x(:), y(:), alpha(:), betamax, xpa, ypa
      integer (I4B) :: minloc
      intent (in) :: x, y, alpha, betamax, minloc
      intent (out) :: xpa, ypa
      end subroutine offcentre
end interface


interface
      subroutine outdxf(xpoint,ypoint,zki,lzrb,laytxt,anzzwg,           &
     &  xgeb,ygeb,txtgeb,txthoc,anzgeb,lplayt)
      use femtypes
      implicit none
      real (DP) :: xpoint(:), ypoint(:), xgeb(:), ygeb(:), txthoc
      integer (I4B) :: lzrb(:), lplayt(:), zki(:,:), anzzwg, anzgeb
      integer (I4B) :: unitid
      character (len=*) :: txtgeb(:), laytxt(:)
      intent (in) :: xpoint, ypoint, zki, lzrb, laytxt, anzzwg
      intent (in) :: xgeb, ygeb, txtgeb, txthoc, anzgeb, lplayt
      end subroutine outdxf
end interface


interface outerprod
      pure function outerprod_c(a,b)
      use femtypes
      implicit none
      complex (DPC) a(:), b(:)
      complex (DPC) outerprod_c(size(a),size(b))
      intent (in) :: a, b
      end function outerprod_c

      pure function outerprod_d(a,b)
      use femtypes
      implicit none
      real (DP) :: a(:), b(:)
      real (DP) :: outerprod_d(size(a),size(b))
      intent (in) :: a, b
      end function outerprod_d
end interface outerprod


interface
      subroutine padapt(ext,int,epsgl)
      use femtypes
      implicit none
      logical, intent (in) :: ext, int
      real (DP), intent (out) :: epsgl
      end subroutine padapt
end interface


interface palloc
      subroutine palloc_dp(ptr,n)
      use femtypes
      implicit none
      integer(I4B) :: n
      real(DP), pointer:: ptr(:)
      intent(in) :: n
      end

      subroutine palloc_sp(ptr,n)
      use femtypes
      implicit none
      integer(I4B) :: n
      real(SP), pointer:: ptr(:)
      intent(in) :: n
      end

      subroutine palloc_dpc(ptr,n)
      use femtypes
      implicit none
      integer(I4B) :: n
      complex(DPC), pointer:: ptr(:)
      intent(in) :: n
      end

      subroutine palloc_spc(ptr,n)
      use femtypes
      implicit none
      integer(I4B) :: n
      complex(SPC), pointer:: ptr(:)
      intent(in) :: n
      end

      subroutine palloc_i(ptr,n)
      use femtypes
      implicit none
      integer(I4B) :: n
      integer(I4B), pointer:: ptr(:)
      intent(in) :: n
      end

      subroutine palloc_l(ptr,n)
      use femtypes
      implicit none
      integer(I4B) :: n
      logical, pointer:: ptr(:)
      intent(in) :: n
      end

      subroutine palloc_ptri(ptr,n)
      use femtypes
      implicit none
      integer(I4B) :: n
      type(ARRPTRI), pointer:: ptr(:)
      intent(in) :: n
      end

      subroutine palloc_dp2(ptr,m,n)
      use femtypes
      implicit none
      integer(I4B) :: m, n
      real(DP), pointer:: ptr(:,:)
      intent(in) :: m, n
      end

      subroutine palloc_sp2(ptr,m,n)
      use femtypes
      implicit none
      integer(I4B) :: m, n
      real(SP), pointer:: ptr(:,:)
      intent(in) :: m, n
      end

      subroutine palloc_dpc2(ptr,m,n)
      use femtypes
      implicit none
      integer(I4B) :: m, n
      complex(DPC), pointer:: ptr(:,:)
      intent(in) :: m, n
      end

      subroutine palloc_spc2(ptr,m,n)
      use femtypes
      implicit none
      integer(I4B) :: m, n
      complex(SPC), pointer:: ptr(:,:)
      intent(in) :: m, n
      end

      subroutine palloc_i2(ptr,m,n)
      use femtypes
      implicit none
      integer(I4B) :: m, n
      integer(I4B), pointer:: ptr(:,:)
      intent(in) :: m, n
      end

      subroutine palloc_l2(ptr,m,n)
      use femtypes
      implicit none
      integer(I4B) :: m, n
      logical, pointer:: ptr(:,:)
      intent(in) :: m, n
      end
end interface


interface
      subroutine pardiso_solver(a,b,x,n,ia,ja)
      use femtypes
      implicit none
      complex (DPC), pointer :: a(:), b(:)
      complex (DPC) x(:)
      integer (I4B), pointer :: ia(:),ja(:)
      integer (I4B) n
      intent (in) :: n
      intent (out) :: x
      end subroutine pardiso_solver
end interface


interface
      subroutine pdecoeff(elem,xs,ys,nu,gamma,alpha,beta,f1,f2,em)
      use femtypes
      implicit none
      integer (I4B) elem
      real (DP) xs, ys
      complex (DPC) :: nu(:,:,:,:), gamma(:,:,:)
      complex (DPC) :: beta(:,:,:), alpha(:,:), f1(:), f2(:,:)
      complex (DPC), optional :: em(:,:)
      intent (in) :: elem, xs, ys
      intent (out) :: nu, gamma, beta, alpha, f1, f2, em
      end subroutine pdecoeff
end interface


interface
      subroutine petsceigen(a, n, eps, ia, ja)
      use femtypes
      implicit none
      complex (DPC), pointer :: a(:)
      integer (I4B)  n 
      integer (I4B), pointer :: ia(:), ja(:)
      real (DP) eps
      intent (in) :: n, eps
      end subroutine petsceigen
end interface

!interface petscsolver
!      subroutine petscsolver_DPC(a, b, x, n, eps, ia, ja)
!      use femtypes
 !     implicit none
  !    complex (DPC), pointer :: a(:)
   !   complex (DPC) b(:), x(:)
  !    integer (I4B)  n 
  !    integer (I4B), pointer :: ia(:), ja(:)
  !    real (DP) eps
  !    intent (in) :: n, eps
  !    intent (inout) :: b, x
  !    end subroutine petscsolver_DPC!
!
!      subroutine petscsolver_DP(a, b, x, n, eps, ia, ja)
!      use femtypes
!      implicit none
!      real (DP), pointer :: a(:)
!      real (DP) b(:), x(:)
!      integer (I4B)  n 
!      integer (I4B), pointer :: ia(:), ja(:)
!      real (DP) eps
!      intent (in) :: n, eps
!      intent (inout) :: b, x
!      end subroutine petscsolver_DP
!end interface

interface
      subroutine petscsolver(a, b, x, n, eps, ia, ja)
      use femtypes
      implicit none
      real (DP), pointer :: a(:)
      real (DP) b(:), x(:)
      integer (I4B)  n 
      integer (I4B), pointer :: ia(:), ja(:)
      real (DP) eps
      intent (in) :: n, eps
      intent (inout) :: b, x
      end subroutine petscsolver
end interface

interface
      subroutine pfeil(x0,y0,lenx,leny,pfwink,pfkopf,head,mitte)
      use femtypes
      implicit none
      real (DP) x0,y0,lenx,leny,pfwink,pfkopf
      logical head,mitte
      intent (in) :: x0,y0,lenx,leny,pfwink,pfkopf,head,mitte
      end subroutine pfeil
end interface

 
interface
      subroutine pfeil2(x0,y0,lenx,leny,pfwink,pfkopf,pfbrt,mitte)
      use femtypes
      implicit none
      real (DP) x0,y0,lenx,leny,pfwink,pfkopf,pfbrt
      logical mitte
      intent (in) :: x0,y0,lenx,leny,pfwink,pfkopf,pfbrt,mitte
      end subroutine pfeil2
end interface


interface
      subroutine physicsinfo(physics,nnat,dimensions)
      use femtypes
      implicit none
      integer (I4B) :: nnat, dimensions
      character (len=*) :: physics
      intent (in) :: physics, dimensions
      intent (out) :: nnat
      end subroutine physicsinfo
end interface


interface
      logical function pinnen(bzi,start,ende,zki,xbk,ybk,x,y)
      use femtypes
      implicit none
      integer (I4B) start,ende,bzi(:),zki(:,:)
      real (DP) xbk(:),ybk(:),x,y
      intent (in) :: bzi,start,ende,zki,xbk,ybk,x,y
      end function pinnen
end interface


interface
      subroutine pktply(j,pzahl,kanz,kknot,xm,ym,xn,yn,neg)
      use femtypes
      implicit none
      integer (I4B) j, pzahl, kanz(:), kknot(:)
      real (DP)  xn(:), yn(:), xm, ym
      logical neg
      intent (in) :: j, pzahl, kanz, kknot, xm, ym, xn, yn
      intent (out) :: neg
      end subroutine pktply
end interface


interface
      subroutine pline(accur,fakt1,scalfk,laytxt,layrb,layanz,xpoint,   &
     &  ypoint,zki,lzrb,zpz,length,anzzwg,anzknt,okall,txtlen)
      use femtypes
      implicit none
      integer (I4B) :: anzzwg, anzknt, txtlen, layanz
      integer (I4B), pointer :: lzrb(:), zki(:,:), zpz(:), layrb(:,:)
      real (DP) :: accur, fakt1, scalfk
      real (DP), pointer :: xpoint(:), ypoint(:), length(:)
      logical :: okall
      character (len=txtlen), pointer :: laytxt(:)
      intent (in) :: accur, fakt1, scalfk, txtlen
      intent (out) :: xpoint, ypoint, zki, lzrb, zpz, length, okall
      intent (inout) :: laytxt, layrb, layanz, anzzwg, anzknt
      end subroutine pline
end interface


interface
      subroutine plottensor(px,py,xa,xe,ya,ye)
      use femtypes
      implicit none
      integer (I4B) px, py
      real (DP) xa, xe, ya, ye
      intent (in) :: px, py, xa, xe, ya, ye
      end subroutine plottensor
end interface


interface
      function point_loss(ielem,lambda)
      use femtypes
      implicit none
      integer (I4B) :: ielem
      real (DP) :: lambda(3), point_loss
      intent (in) :: ielem, lambda
      end function point_loss
end interface


interface
      subroutine posarc(xbk,ybk,zki,neukzi,xneu,yneu,x1,y1,x2,y2,      &
     &  x3,y3,x4,y4,neighb,err)
      use femtypes
      implicit none
      integer (I4B) zki(:,:),neukzi, err
      real (DP) xbk(:),ybk(:),xneu,yneu,x1,y1,x2,y2, x3, y3, x4, y4
      logical neighb
      intent (in) :: xbk, ybk, zki,neukzi, x1, y1, x2, y2
      intent (in) :: x3, y3, x4, y4, neighb
      intent (out) :: xneu, yneu, err
      end subroutine posarc
end interface


interface
      subroutine poyntg(x,xn,yn,e,en,n,matzif,omega,gbz,geb,reell)
      use femtypes
      implicit none
      real (DP) xn(:), yn(:)
      complex (DPC) x(:)
      real (DP) omega
      integer (I4B) e(:,:), en(:,:), n, matzif(:), geb(:), gbz
      logical reell
      intent (in) :: x, xn, yn, e, en, n, matzif, omega, gbz, geb, reell
      end subroutine poyntg
end interface


interface
      subroutine preassemb(oldep)
      use femtypes
      implicit none
      integer (I4B), optional, pointer:: oldep(:,:)
      end subroutine preassemb
end interface


interface
      subroutine principlevalues(vx,vy,vxy,lambda1,lambda2,angle)
      use femtypes
      implicit none
      real (DP) :: vx, vy, vxy, lambda1, lambda2
      real (DP), optional :: angle
      intent (in) :: vx, vy, vxy
      intent (out) :: lambda1, lambda2, angle
      end subroutine principlevalues
end interface

interface
      subroutine print_error(error_level,msg,warning)
      use femtypes
      implicit none
      integer     (I4B)               :: error_level
      character (len=*),     optional :: warning
      character (len=*)               :: msg
      end subroutine print_error
end interface


interface
      subroutine putlin(x1,y1,x2,y2,r,phi1,phi2,gerad,accur,fakt1,      &
     &  layer,ok,xpoint,ypoint,zki,lzrb,zpz,length,anzzwg,anzknt)
      use femtypes
      implicit none
      integer (I4B) :: layer, anzzwg, anzknt
      integer (I4B), pointer :: lzrb(:), zki(:,:), zpz(:)
      real (DP) :: x1, x2, y1, y2, r, phi1, phi2, accur, fakt1
      real (DP), pointer :: xpoint(:), ypoint(:), length(:)
      logical :: gerad, ok
      intent (in) :: x1, y1, x2, y2, r, phi1, phi2, gerad, layer, accur
      intent (in) :: fakt1
      intent (out) :: length, ok
      intent (inout) :: xpoint, ypoint, zki, lzrb, zpz, anzzwg, anzknt
      end subroutine putlin
end interface


interface
      subroutine putmatparam(iparam,str)
      use femtypes
      implicit none
      integer (I4B) iparam
      character (len=*) str
      intent (in) :: iparam,str
      end subroutine putmatparam
end interface


interface
      subroutine putnod(x,y,accur,node,xpoint,ypoint,anzknt)
      use femtypes
      implicit none
      real (DP) x,y,accur,xpoint(:),ypoint(:)
      integer (I4B) node,anzknt
      intent (in) :: x, y, accur
      intent (out) :: node
      intent (inout) :: xpoint, ypoint, anzknt
      end subroutine putnod
end interface


interface putpostsetting
      subroutine putpostsetting_c(key,value)
      use femtypes
      implicit none
      character (len=*) key
      character (len=*) value
      intent (in) key
      intent (in) value
      end subroutine putpostsetting_c

      subroutine putpostsetting_f(key,value)
      use femtypes
      implicit none
      character (len=*) key
      real (DP) value
      intent (in) key
      intent (in) value
      end subroutine putpostsetting_f

      subroutine putpostsetting_i(key,value)
      use femtypes
      implicit none
      character (len=*) key
      integer (I4B) value
      intent (in) key
      intent (in) value
      end subroutine putpostsetting_i
end interface


interface putsetting
      subroutine putsetting_c(key,value)
      use femtypes
      implicit none
      character (len=*) key
      character (len=*) value
      intent (in) key
      intent (in) value
      end subroutine putsetting_c

      subroutine putsetting_i(key,value)
      use femtypes
      implicit none
      character (len=*) key
      real (DP) value
      intent (in) key
      intent (in) value
      end subroutine putsetting_i

      subroutine putsetting_f(key,value)
      use femtypes
      implicit none
      character (len=*) key
      integer (I4B) value
      intent (in) key
      intent (in) value
      end subroutine putsetting_f
end interface


interface qmrsolver
      subroutine qmrsolver_DPC(lower,upper,diag,b,x,n,eps,ia,ja,epsgl,resgl,symm)
      use femtypes
      implicit none
      complex (DPC) lower(:),upper(:),diag(:),b(:),x(:)
      integer (I4B) n,ia(:),ja(:)
      real (DP) eps,epsgl,resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: lower, upper, diag, b, x
      end subroutine qmrsolver_DPC

      subroutine qmrsolver_DP(lower,upper,diag,b,x,n,eps,ia,ja,epsgl,resgl,symm)
      use femtypes
      implicit none
      real (DP) lower(:),upper(:),diag(:),b(:),x(:)
      integer (I4B) n,ia(:),ja(:)
      real (DP) eps,epsgl,resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: lower, upper, diag, b, x
      end subroutine qmrsolver_DP
end interface


interface qsort
      subroutine qsort_dp(arr,n)
      use femtypes
      implicit none
      real (DP) arr(:)
      integer (I4B) n
      intent (in) :: n
      intent (inout) :: arr
      end subroutine qsort_dp
end interface


interface qsortindex
      subroutine qsortindex_dp(arr,indx,n)
      use femtypes
      implicit none
      real (DP) :: arr(:)
      integer (I4B) :: indx(:), n
      intent (in) :: arr, n
      intent (out) :: indx
      end subroutine qsortindex_dp
      
      subroutine qsortindex_sp(arr,indx,n)
      use femtypes
      implicit none
      real (SP) :: arr(:)
      integer (I4B) :: indx(:), n
      intent (in) :: arr, n
      intent (out) :: indx
      end subroutine qsortindex_sp

      subroutine qsortindex_i4(arr,indx,n)
      use femtypes
      implicit none
      integer (I4B) :: arr(:)
      integer (I4B) :: indx(:), n
      intent (in) :: arr, n
      intent (out) :: indx
      end subroutine qsortindex_i4
end interface


interface
      subroutine rbtxt(txtgeb,anzgeb,alrb,btrb,donetx,lplayt,layrb,nnat)
      use femtypes
      implicit none
      integer (I4B) :: anzgeb, lplayt(:), nnat
      integer (I4B) :: layrb(:,:), donetx(:)
      complex (DPC) :: alrb(:,:), btrb(:,:)
      character (len=*) :: txtgeb(:)
      intent (in) :: txtgeb, anzgeb, lplayt, nnat
      intent (out) ::  alrb, btrb, donetx, layrb
      end subroutine rbtxt
end interface


interface
      subroutine rcm(root, xadj, adjncy, mask, perm, ccsize, offset, cm)
      use femtypes
      implicit none
      integer (I4B) root, xadj(:), adjncy(:), mask(:), perm(:)
      integer (I4B) ccsize, offset
      logical cm
      intent (in) :: root, adjncy, offset, cm
      intent (out) :: perm, ccsize
      intent (inout) :: xadj, mask
      end subroutine rcm
end interface

interface
      subroutine read_femsettingsfile(path)
      use femtypes
      implicit none
      character (len=200)  :: path
      end subroutine
end interface

interface
      subroutine readky(key,string,float,int,line)
      use femtypes
      implicit none
      real (DP) :: float
      integer (I4B) :: line, int, key
      character (len=*) :: string
      intent (out) :: key, string, float, int, line
      end subroutine readky
end interface


interface
      subroutine readmat(matname,found,matindex)
      use femtypes
      implicit none
      character (len=*) matname
      integer (I4B) matindex
      logical found
      intent (out) :: found, matindex
      intent (in) :: matname
      end subroutine readmat
end interface


interface
      subroutine readnetin(gbz,gzz,gkz,zki,kzrb,zrb,xbk,ybk,matzif,ok,        &
     &                     alrb,btrb,nnat)
      use femtypes
      implicit none
      integer (I4B) gbz, gzz, gkz, nnat
      integer (I4B), pointer ::  kzrb(:,:), zrb(:,:), matzif(:),  zki(:,:)
      real (DP), pointer ::  xbk(:), ybk(:)
      complex (DPC), pointer :: alrb(:,:), btrb(:,:)
      logical ok
      intent (in ) :: nnat
      intent (out) :: gbz, gzz, gkz, zki, kzrb, zrb, xbk, ybk, matzif
      intent (out) :: ok, alrb, btrb
      end subroutine readnetin
end interface


interface
      subroutine readphysics(physics, dimensions, nnat, ok)
      use femtypes
      implicit none
      integer (I4B) :: dimensions, nnat
      character (len=*) :: physics
      logical :: ok
      intent (in) :: physics, dimensions
      intent (out) :: nnat, ok
      end subroutine readphysics
end interface


interface
      subroutine readpostsetting(callback)
      use femtypes
      implicit none
      character (len=*):: callback
      end subroutine readpostsetting
end interface


interface
      subroutine readsetting
      use femtypes
      implicit none
      end subroutine readsetting
end interface


interface realloc
      subroutine realloc_dp(array, n)
      use femtypes
      implicit none
      integer(I4B) :: n
      real(DP), allocatable :: array(:)
      intent(in) :: n
      intent(inout) :: array
      end subroutine realloc_dp

      subroutine realloc_dp_inout(arrayin, n, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n
      real(DP), allocatable :: arrayin(:)
      real(DP), allocatable :: arrayout(:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
      end subroutine realloc_dp_inout

      subroutine realloc_dp2(array, n, m)
      use femtypes
      implicit none
      integer(I4B) :: n, m
      real(DP), allocatable :: array(:,:)
      intent(in) :: n
      intent(inout) :: array
      end subroutine realloc_dp2

      subroutine realloc_dp2_inout(arrayin, n, m, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n, m
      real(DP), allocatable :: arrayin(:,:)
      real(DP), allocatable :: arrayout(:,:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
      end subroutine realloc_dp2_inout

      subroutine realloc_sp(array, n)
      use femtypes
      implicit none
      integer(I4B) :: n
      real(SP), allocatable :: array(:)
      intent(in) :: n
      intent(inout) :: array
      end subroutine realloc_sp

      subroutine realloc_sp_inout(arrayin, n, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n
      real(SP), allocatable :: arrayin(:)
      real(SP), allocatable :: arrayout(:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
      end subroutine realloc_sp_inout

      subroutine realloc_dpc(array, n)
      use femtypes
      implicit none
      integer(I4B) :: n
      complex(DPC), allocatable :: array(:)
      intent(in) :: n
      intent(inout) :: array
      end subroutine realloc_dpc

      subroutine realloc_dpc_inout(arrayin, n, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n
      complex(DPC), allocatable :: arrayin(:)
      complex(DPC), allocatable :: arrayout(:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
     end subroutine realloc_dpc_inout

      subroutine realloc_spc(array, n)
      use femtypes
      implicit none
      integer(I4B) :: n
      complex(SPC), allocatable :: array(:)
      intent(in) :: n
      intent(inout) :: array
      end subroutine realloc_spc

      subroutine realloc_spc_inout(arrayin, n, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n
      complex(SPC), allocatable :: arrayin(:)
      complex(SPC), allocatable :: arrayout(:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
      end subroutine realloc_spc_inout

      subroutine realloc_i(array, n)
      use femtypes
      implicit none
      integer(I4B) :: n
      integer(I4B), allocatable :: array(:)
      intent(in) :: n
      intent(inout) :: array
      end subroutine realloc_i

      subroutine realloc_i_inout(arrayin, n, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n
      integer(I4B), allocatable :: arrayin(:)
      integer(I4B), allocatable :: arrayout(:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
      end subroutine realloc_i_inout

      subroutine realloc_i2(array, n, m)
      use femtypes
      implicit none
      integer(I4B) :: n, m
      integer(I4B), allocatable :: array(:,:)
      intent(in) :: n
      intent(inout) :: array
      end subroutine realloc_i2

      subroutine realloc_i2_inout(arrayin, n, m, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n, m
      integer(I4B), allocatable :: arrayin(:,:)
      integer(I4B), allocatable :: arrayout(:,:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
      end subroutine realloc_i2_inout

      subroutine realloc_l(array, n)
      use femtypes
      implicit none
      integer(I4B) :: n
      logical, allocatable :: array(:)
      intent(in) :: n
      intent(inout) :: array
      end subroutine realloc_l

      subroutine realloc_l_inout(arrayin, n, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n
      logical, allocatable :: arrayin(:)
      logical, allocatable :: arrayout(:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
      end subroutine realloc_l_inout

      subroutine realloc_l2(array, n, m)
      use femtypes
      implicit none
      integer(I4B) :: n, m
      logical, allocatable :: array(:,:)
      intent(in) :: n
      intent(inout) :: array
      end subroutine realloc_l2

      subroutine realloc_l2_inout(arrayin, n, m, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n, m
      logical, allocatable :: arrayin(:,:)
      logical, allocatable :: arrayout(:,:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
      end subroutine realloc_l2_inout

      subroutine realloc_ch(array, n, txtlen)
      use femtypes
      implicit none
      integer(I4B) :: n, txtlen
      character (len=txtlen), allocatable :: array(:)
      intent(in) :: n, txtlen
      intent(inout) :: array
      end subroutine realloc_ch

      subroutine realloc_ch_inout(arrayin, n, txtlen, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n, txtlen
      character (len=txtlen), allocatable :: arrayin(:)
      character (len=txtlen), allocatable :: arrayout(:)
      intent(in) :: n, txtlen
      intent(out) :: arrayout
      intent(inout) :: arrayin
      end subroutine realloc_ch_inout
end interface


interface reallocate
      function reallocate_dp(p,n)
      use femtypes
      implicit none
      real (DP), pointer :: p(:), reallocate_dp(:)
      integer (I4B) n
      intent (in) :: n
      end function reallocate_dp

      function reallocate_dp2(p,n,m)
      use femtypes
      implicit none
      real (DP), pointer :: p(:,:), reallocate_dp2(:,:)
      integer (I4B)      :: n, m
      intent(in)         :: n, m
      end function reallocate_dp2

      function reallocate_sp(p,n)
      use femtypes
      implicit none
      real (SP), pointer :: p(:), reallocate_sp(:)
      integer (I4B) n
      intent (in) :: n
      end function reallocate_sp

      function reallocate_dpc(p,n)
      use femtypes
      implicit none
      complex (DPC), pointer :: p(:), reallocate_dpc(:)
      integer (I4B) n
      intent (in) :: n
      end function reallocate_dpc

      function reallocate_spc(p,n)
      use femtypes
      implicit none
      complex (SPC), pointer :: p(:), reallocate_spc(:)
      integer (I4B) n
      intent (in) :: n
      end function reallocate_spc

      function reallocate_i(p,n)
      use femtypes
      implicit none
      integer (I4B), pointer :: p(:), reallocate_i(:)
      integer (I4B) n
      intent (in) :: n
      end function reallocate_i

      function reallocate_i2(p,n,m)
      use femtypes
      implicit none
      integer (I4B), pointer :: p(:,:), reallocate_i2(:,:)
      integer (I4B) n, m
      intent (in) :: n, m
      end function reallocate_i2

      function reallocate_ch(p,n,txtlen)
      use femtypes
      implicit none
      integer (I4B) n, txtlen
      character (len=txtlen), pointer :: p(:), reallocate_ch(:)
      intent (in) :: n, txtlen
      end function reallocate_ch

      function reallocate_l(p,n)
      use femtypes
      implicit none
      logical, pointer :: p(:), reallocate_l(:)
      integer (I4B) n
      intent(in) :: n
      end function reallocate_l

      function reallocate_l2(p,n,m)
      use femtypes
      implicit none
      logical, pointer :: p(:,:), reallocate_l2(:,:)
      integer (I4B) n, m
      intent(in) :: n, m
      end function reallocate_l2
      
      function reallocate_perr(p,n)
      use femtypes
      implicit none
      type(ARRPTRI), pointer :: p(:), reallocate_perr(:)
      integer (I4B) n
      intent(in) :: n
      end function reallocate_perr
end interface


interface
      subroutine redmat(matname,imat,found)
      use femtypes
      implicit none
      integer (I4B) imat
      character(len=*) matname 
      logical found
      intent (in) :: matname
      intent (out) :: found, imat
      end subroutine redmat
end interface


interface
      subroutine reduba(e,en,n,p,xn,yn,x,kzi)
      use femtypes
      implicit none
      integer (I4B) e(:,:),en(:,:),n,p
      integer (I4B) kzi(:)
      real (DP) xn(:),yn(:)
      complex (DPC) x(:)
      intent (in) :: en, n, p
      intent (inout) :: e, xn, yn, x, kzi
      end subroutine reduba
end interface


interface
      subroutine reduce(e,en,n,p,xn,yn,x,kzi)
      use femtypes
      implicit none
      integer (I4B)  e(:,:), en(:,:), n, p, kzi(:)
      real (DP) xn(:),yn(:)
      complex (DPC) x(:)
      intent (in) :: en, n, p
      intent (inout) e, xn, yn, x, kzi
      end subroutine reduce
end interface


interface
      subroutine regionlist(bzi,bzip,bzil,zpz,branchpt,brnodes,region,  &
     &  head,nbds,numnodes)
      use femtypes
      use triangtypes
      implicit none
      integer (I4B) bzi(:), bzip(:), bzil(:), zpz(:)
      integer (I4B) branchpt(:), brnodes(:), region, nbds
      integer (I4B), pointer :: numnodes(:)
      type (arrpt), pointer :: head(:)
      intent (in) :: bzi, bzip, bzil, zpz, branchpt, brnodes, region
      intent (out) :: nbds
      end subroutine regionlist
end interface


interface
      subroutine resid_unscaled(lower,upper,diag,b,x,ia,ja,n,res)
      use femtypes
      implicit none
      complex (DPC) lower(:), upper(:), diag(:), x(:), b(:), res(:)
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: upper, lower, diag, b, x, ia, ja, n
      intent (out) :: res
      end subroutine resid_unscaled
end interface


interface resid
      subroutine resid_c(lower,upper,diag,b,x,ia,ja,n,res)
      use femtypes
      implicit none
      complex (SPC) lower(:),upper(:),diag(:),x(:),b(:),res(:)
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: upper, lower, diag, b, x, ia, ja, n
      intent (out) :: res
      end subroutine resid_c

      subroutine resid_z(lower,upper,diag,b,x,ia,ja,n,res)
      use femtypes
      implicit none
      complex (DPC) lower(:),upper(:),diag(:),x(:),b(:),res(:)
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: upper, lower, diag, b, x, ia, ja, n
      intent (out) :: res
      end subroutine resid_z

      subroutine resid_d(lower,upper,diag,b,x,ia,ja,n,res)
      use femtypes
      implicit none
      real (DP) lower(:),upper(:),diag(:),x(:),b(:),res(:)
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: upper, lower, diag, b, x, ia, ja, n
      intent (out) :: res
      end subroutine resid_d

      subroutine resid_s(lower,upper,diag,b,x,ia,ja,n,res)
      use femtypes
      implicit none
      real (SP) lower(:),upper(:),diag(:),x(:),b(:),res(:)
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: upper, lower, diag, b, x, ia, ja, n
      intent (out) :: res
      end subroutine resid_s
end interface


interface resid_unscaled
      subroutine resid_unscaled_c(lower,upper,diag,b,x,ia,ja,n,res)
      use femtypes
      implicit none
      complex (SPC) lower(:), upper(:), diag(:), x(:), b(:), res(:)
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: upper, lower, diag, b, x, ia, ja, n
      intent (out) :: res
      end subroutine resid_unscaled_c

      subroutine resid_unscaled_z(lower,upper,diag,b,x,ia,ja,n,res)
      use femtypes
      implicit none
      complex (DPC) lower(:), upper(:), diag(:), x(:), b(:), res(:)
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: upper, lower, diag, b, x, ia, ja, n
      intent (out) :: res
      end subroutine resid_unscaled_z
end interface


interface residual
      subroutine residual(ext,int,matvar,errcode,res,sumres,sumref)
      use femtypes
      implicit none
      integer (I4B) errcode
      real (DP) :: res(:,:), sumres(:), sumref(:)
      logical :: ext, int, matvar
      intent (in) :: ext, int, matvar
      intent (out) :: errcode, res, sumres, sumref
      end subroutine residual
end interface


interface
      subroutine revis()
      use femtypes
      implicit none
      end subroutine revis
end interface


interface
      subroutine rgauss(ra,rb,n,rx)
      use femtypes
      implicit none
      integer (I4B) n
      real (DP) ra(:,:),rb(:),rx(:)
      intent (in) :: ra,rb,n
      intent (out) :: rx
      end subroutine rgauss
end interface


interface
      subroutine rootls(root, xadj, adjncy, mask, nlvl, xls, ls, offset)
      use femtypes
      implicit none
      integer (I4B) root, xadj(:), adjncy(:), mask(:), nlvl
      integer (I4B) xls(:), ls(:), offset
      intent (in) :: root, xadj, adjncy, offset
      intent (out) :: nlvl, xls, ls
      intent (inout) :: mask
      end subroutine rootls
end interface


interface
      subroutine rotatehom(a,theta,axis)
      use femtypes
      implicit none
      real (DP) a(4,4), theta
      integer (I4B) axis
      intent (in) theta, axis
      intent (inout) a
      end subroutine rotatehom
end interface


interface scal
      subroutine scal_c(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
      use femtypes
      implicit none
      complex (SPC) upper(:), lower(:), diag(:), x(:), b(:)
      real (SP) bscal
      logical symm
      integer (I4B) ia(:), ja(:), n
      intent (in) :: ia, ja, n, symm
      intent (out) :: bscal
      intent (inout) :: lower, upper, diag, b, x
      end subroutine scal_c

      subroutine scal_z(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
      use femtypes
      implicit none
      complex (DPC) upper(:), lower(:), diag(:), x(:), b(:)
      real (DP) bscal
      integer (I4B) ia(:), ja(:), n
      logical symm
      intent (in) :: ia, ja, n, symm
      intent (out) :: bscal
      intent (inout) :: lower, upper, diag, b, x
      end subroutine scal_z

      subroutine scal_d(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
      use femtypes
      implicit none
      real (DP) upper(:), lower(:), diag(:), x(:), b(:), bscal
      integer (I4B) ia(:), ja(:), n
      logical symm
      intent (in) :: ia, ja, n, symm
      intent (out) :: bscal
      intent (inout) :: lower, upper, diag, b, x
      end subroutine scal_d

      subroutine scal_s(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
      use femtypes
      implicit none
      real (SP) upper(:), lower(:), diag(:), x(:), b(:), bscal
      integer (I4B) ia(:), ja(:), n
      logical symm
      intent (in) :: ia, ja, n, symm
      intent (out) :: bscal
      intent (inout) :: lower, upper, diag, b, x
      end subroutine scal_s
end interface


interface
      subroutine scalehom(a,sx,sy,sz)
      use femtypes
      implicit none
      real (DP) a(4,4), sx, sy, sz
      intent (in) sx, sy, sz
      intent (inout) a
      end subroutine scalehom
end interface


interface
      subroutine scatterparameters(xa,ya,xe,ye,swr,rc,tc)
      use femtypes
      implicit none
      real (DP) :: xa, xe, ya, ye
      real (DP) :: swr, rc, tc
      intent (in) :: xa, ya, xe, ye
      intent (out) :: swr, rc, tc
      end subroutine scatterparameters
end interface


interface
      logical function schn(x1,y1,x2,y2,x3,y3,x4,y4)
      use femtypes
      implicit none
      real (DP) x1,y1,x2,y2,x3,y3,x4,y4
      intent (in) :: x1,y1,x2,y2,x3,y3,x4,y4
      end function schn
end interface


interface
      subroutine schn1(x1,y1,x2,y2,x3,y3,x4,y4,xs,ys,schni)
      use femtypes
      implicit none
      real (DP) x3,y3,x4,y4,x1,y1,x2,y2,xs,ys
      logical schni
      intent (in) :: x1, y1, x2, y2, x3, y3, x4, y4
      intent (out) :: xs, ys, schni
      end subroutine schn1
end interface


interface
      subroutine schnit(x1,y1,x2,y2,x3,y3,x4,y4,t1,t2,ok)
      use femtypes
      implicit none
      real (DP) x1, y1, x2, y2, x3, y3, x4, y4, t1, t2
      logical ok
      intent (in) :: x1, y1, x2, y2, x3, y3, x4, y4
      intent (out) :: t1, t2, ok
      end subroutine schnit
end interface


interface
      subroutine second(ttime)
      use femtypes
      implicit none
      real (DP) ttime
      intent (out) :: ttime      
      end subroutine second
end interface


interface
      subroutine sectio(sction)
      use femtypes
      implicit none
      integer (I4B) sction
      intent (out) :: sction
      end subroutine sectio
end interface


interface
      subroutine setprocesspriority(priority)
      use femtypes
      implicit none
      integer (I4B) priority
      intent (in) :: priority
      end subroutine setprocesspriority
end interface


interface
      subroutine setstandardvalues
      use femtypes
      implicit none
      end subroutine setstandardvalues
end interface


interface
      subroutine setze(ellist,zahle,zahlk,elem)
      use femtypes
      implicit none
      integer (I4B) ellist(:), zahle, zahlk, elem
      intent (in) :: elem
      intent (inout) :: ellist,zahle,zahlk
      end subroutine setze
end interface


interface
      subroutine setzeg(ellist,zahle,zahlk,elem,node)
      use femtypes
      implicit none
      integer (I4B) ellist(:), zahle, zahlk, elem, node
      intent (in) :: elem, node
      intent (inout) :: ellist, zahle, zahlk
      end subroutine setzeg
end interface


interface
      subroutine setzwg(length,anzzwg,anzknt,zki,zpz,zpp)
      use femtypes
      implicit none
      integer (I4B) anzknt,anzzwg,zki(:,:),zpz(:)
      real (DP) length(:),zpp(:)
      intent (in) :: length,anzzwg,anzknt,zki
      intent (out) :: zpz,zpp
      end subroutine setzwg
end interface


interface
      pure subroutine shape(l, polylo, polyhi, vec)
      use femtypes
      implicit none
      integer (I4B) polylo, polyhi
      real (DP) l(:)
      real (DP) vec(:)
      intent (in) :: l, polylo, polyhi
      intent (out) :: vec
      end subroutine shape
end interface


interface      
      subroutine shapea(l, polylo, polyhi, vec)
      use femtypes
      implicit none
      integer (I4B) polylo, polyhi
      real (DP) l(3)
      real (DP) vec(:)
      intent (in) :: l, polylo, polyhi
      intent (out) :: vec
      end subroutine shapea
end interface


interface
      pure subroutine shapefunction(lambda,xl,yl,polylo,polyhi,nff,          &
     &  gradient,xsi,gxsi,errcode)
      use femtypes
      implicit none
      integer (I4B) polylo, polyhi, nff, errcode
      real (DP) xl(:), yl(:), lambda(:)
      real (DP) xsi(:), gxsi(:,:)
      logical gradient
      intent (in) :: lambda, xl, yl, polylo, polyhi, nff, gradient
      intent (out) :: xsi, gxsi, errcode
      end subroutine shapefunction
end interface


interface      
      subroutine shapew(l, polylo, polyhi, vec)
      use femtypes
      implicit none
      integer (I4B) polylo, polyhi
      real (DP) l(3)
      real (DP) vec(:)
      intent (in) :: l, polylo, polyhi
      intent (out) :: vec
      end subroutine shapew
end interface


interface
      subroutine snkrgr(xa,ya,xb,yb,xc,yc,xd,yd,xm,ym,rl1,rl2,t1,t2,err)
      use femtypes
      implicit none
      real (DP) xa,ya,xb,yb,xc,yc,xd,yd,xm,ym,rl1,rl2,t1,t2
      integer (I4B) err
      intent (in) :: xa, ya, xb, yb, xc, yc, xd, yd, xm, ym
      intent (out) :: rl1, rl2, t1, t2, err
      end subroutine snkrgr
end interface


interface
      subroutine snkrkr(xa,ya,xb,yb,xc,yc,xd,yd,xm1,ym1,xm2,ym2,        &
     &  t1u,t1v,t2u,t2v,xpu,ypu,xpv,ypv,err)
      use femtypes
      implicit none
      real (DP) xa, ya, xb, yb, xc, yc, xd, yd, xm1, ym1, xm2, ym2
      real (DP) t1u, t1v, t2u, t2v, xpu, ypu, xpv, ypv
      integer (I4B) err
      intent (in) :: xa, ya, xb, yb, xc, yc, xd, yd, xm1, ym1, xm2, ym2
      intent (out) :: t1u, t1v, t2u, t2v, xpu, ypu, xpv, ypv
      end subroutine snkrkr
end interface


interface
      subroutine snzwg(xh,yh,zki,i,j,xp,yp,ls)
      use femtypes
      implicit none
      integer (I4B) zki(:,:),i,j
      real (DP) xh(:),yh(:),xp,yp
      logical ls
      intent (in) :: xh, yh, zki,i, j
      intent (out) :: xp, yp, ls
      end subroutine snzwg
end interface


interface solve2
      subroutine solve2_c(lower,upper,diag,b,x,n,eps,ia,ja,epsgl,resgl,symm)
      use femtypes
      implicit none
      complex (SPC) lower(:), upper(:), diag(:), b(:), x(:)
      integer (I4B) n, ia(:), ja(:)
      real (SP) eps, epsgl, resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: lower, upper, diag, b, x
      end subroutine solve2_c

      subroutine solve2_z(lower,upper,diag,b,x,n,eps,ia,ja,epsgl,resgl,symm)
      use femtypes
      implicit none
      complex (DPC) lower(:), upper(:), diag(:), b(:), x(:)
      integer (I4B) n, ia(:), ja(:)
      real (DP) eps, epsgl, resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: lower, upper, diag, b, x
      end subroutine solve2_z

      subroutine solve2_d(lower,upper,diag,b,x,n,eps,ia,ja,epsgl,resgl,symm)
      use femtypes
      implicit none
      real (DP) lower(:), upper(:), diag(:), b(:), x(:)
      integer (I4B) n, ia(:), ja(:)
      real (DP) eps, epsgl, resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: lower, upper, diag, b, x
      end subroutine solve2_d

      subroutine solve2_s(lower,upper,diag,b,x,n,eps,ia,ja,epsgl,resgl,symm)
      use femtypes
      implicit none
      real (SP) lower(:), upper(:), diag(:), b(:), x(:)
      integer (I4B) n, ia(:), ja(:)
      real (SP) eps, epsgl, resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: lower, upper, diag, b, x
      end subroutine solve2_s
end interface


interface sort_asc_order
      subroutine sort_asc_order_DP(a, ia, ja)
      use femtypes
      implicit none
      real (DP), pointer:: a(:)
      integer (I4B), pointer :: ia(:), ja(:)
      intent(in) :: ia
      intent (inout) :: a, ja
      end subroutine sort_asc_order_DP

      subroutine sort_asc_order_DPC(a, ia, ja)
      use femtypes
      implicit none
      complex (DPC), pointer:: a(:)
      integer (I4B), pointer :: ia(:), ja(:)
      intent(in) :: ia
      intent (inout) :: a, ja
      end subroutine sort_asc_order_DPC
end interface


interface
      subroutine sortnt(anzknt,anzzwg,xpoint,ypoint,zki,zpz,lzrb)
      use femtypes
      implicit none
      integer (I4B) :: anzknt, anzzwg, zki(:,:), zpz(:), lzrb(:)
      real (DP) :: xpoint(:), ypoint(:)
      intent (in) :: anzknt, anzzwg
      intent (inout) :: xpoint, ypoint, zki, zpz, lzrb
      end subroutine sortnt
end interface


interface
      subroutine sorttr(gbz,gzz,bzi,bzip,bzil,zki,zrb,ok)
      use femtypes
      implicit none
      integer (I4B) :: gbz,gzz
      integer (I4B) :: bzi(:), bzip(:), bzil(:), zki(:,:), zrb(:,:)
      logical :: ok
      intent (in) :: gbz, gzz, bzip, zki, zrb
      intent (out) :: ok
      end subroutine sorttr
end interface


interface
      subroutine splitencrelement(i,j,slist,splitted)
      use femtypes
      implicit none
      integer (I4B) i, j, slist(:)
      logical splitted
      intent (in) :: i ,j
      intent (inout) :: slist
      intent (out) splitted
      end subroutine splitencrelement
end interface


interface
      subroutine sptransformpoly(xl,yl,zl,a,n,xo,yo,zo)
      use femtypes
      implicit none
      real (DP) a(4,4)
      real (SP) xl(:), yl(:), zl(:), xo(:), yo(:), zo(:)
      integer (I4B) n
      intent (in) xl, yl, zl, a, n
      intent (out) xo, yo, zo
      end subroutine sptransformpoly
end interface


interface string2number
      subroutine string2number_SP(string,number,ok)
      use femtypes
      implicit none
      character(len=*) :: string
      real (SP) :: number
      logical, optional :: ok
      intent (in) :: string
      intent (out) :: number, ok
      end subroutine string2number_SP

      subroutine string2number_DP(string,number,ok)
      use femtypes
      implicit none
      character(len=*) :: string
      real (DP) :: number
      logical, optional :: ok
      intent (in) :: string
      intent (out) :: number, ok
      end subroutine string2number_DP

      subroutine string2number_I(string,number,ok)
      use femtypes
      implicit none
      character(len=*) :: string
      integer (I4B) :: number
      logical, optional :: ok
      intent (in) :: string
      intent (out) :: number, ok
      end subroutine string2number_I
end interface


interface number2string
      subroutine number2string_SP(number,string,ok)
      use femtypes
      implicit none
      character(len=*)  :: string
      real (SP)         :: number
      logical, optional :: ok
      intent (in)       :: number
      intent (out)      :: string, ok
      end subroutine number2string_SP
      
      subroutine number2string_DP(number,string,ok)
      use femtypes
      implicit none
      character(len=*)  :: string
      real (DP)         :: number
      logical, optional :: ok
      intent (in)       :: number
      intent (out)      :: string, ok
      end subroutine number2string_DP
      
      subroutine number2string_I(number,string,ok)
      use femtypes
      implicit none
      character(len=*)  :: string
      integer (I4B)     :: number
      logical, optional :: ok
      intent (in)       :: number
      intent (out)      :: string, ok
      end subroutine number2string_I
end interface


interface
      subroutine strtok(srce, delim, token)
      use femtypes
      implicit none
      character (len=*) :: token, srce, delim
      intent (in) :: srce, delim
      intent (out) :: token
      end subroutine strtok
end interface


interface
      subroutine strtok1(srce, delim, token)
      use femtypes
      implicit none
      character (len=*) :: token, srce, delim
      intent (in) :: srce, delim
      intent (out) :: token
      end subroutine strtok1
end interface

interface
      subroutine struc(gzz,zki,zrb,xbk,ybk,r,g,b)
      use femtypes
      integer (I4B) gzz,zki(:,:),zrb(:,:)
      real (DP) xbk(:),ybk(:)
      real (SP) r, g, b
      intent (in):: gzz,zki,zrb,xbk,ybk,r,g,b
      end subroutine struc
end interface


interface swap
      pure subroutine swap_c(a,b)
      use femtypes
      implicit none
      complex (DPC) a(:), b(:)
      intent (inout) :: a, b
      end subroutine swap_c

      pure subroutine swap_d(a,b)
      use femtypes
      implicit none
      real (DP) :: a(:), b(:)
      intent (inout) :: a, b
      end subroutine swap_d
end interface swap


interface
      subroutine tab2blank(string)
      use femtypes
      implicit none
      character(len=*) :: string
      intent (inout) :: string
      end subroutine tab2blank
end interface


interface
      subroutine tables(lytabl,layanz,laytxt,layrb,txtlen)
      use femtypes
      implicit none
      integer (I4B) :: layanz, txtlen
      integer (I4B), pointer :: layrb(:,:)
      character (len=txtlen), pointer :: laytxt(:)
      logical :: lytabl
      intent (in) :: txtlen
      intent (out) :: layanz, laytxt, layrb, lytabl
      end subroutine tables
end interface

! TODO: Replace old tensor routines by ne ones
interface ten2mat_r2
      function ten2mat_r2_dpc(ten)
      use femtypes
      implicit none
      complex (DPC) :: ten(3,3), ten2mat_r2_dpc(6)
      intent (in) :: ten
      end function ten2mat_r2_dpc

      function ten2mat_r2_dp(ten)
      use femtypes
      implicit none
      real (DP) :: ten(3,3), ten2mat_r2_dp(6)
      intent (in) :: ten
      end function ten2mat_r2_dp
end interface


interface ten2mat_r3
      function ten2mat_r3_dpc(ten)
      use femtypes
      implicit none
      complex (DPC) :: ten(3,3,3), ten2mat_r3_dpc(6,3)
      intent (in) :: ten
      end function ten2mat_r3_dpc

      function ten2mat_r3_dp(ten)
      use femtypes
      implicit none
      real (DP) :: ten(3,3,3), ten2mat_r3_dp(6,3)
      intent (in) :: ten
      end function ten2mat_r3_dp
end interface


interface ten2mat_r4
      function ten2mat_r4_dpc_old(ten)
      use femtypes
      implicit none
      complex (DPC) :: ten(3,3,3,3), ten2mat_r4_dpc_old(6,6)
      intent (in) :: ten
      end function ten2mat_r4_dpc_old

      function ten2mat_r4_dp_old(ten)
      use femtypes
      implicit none
      real (DP) :: ten(3,3,3,3), ten2mat_r4_dp_old(6,6)
      intent (in) :: ten
      end function ten2mat_r4_dp_old
end interface


interface ten2mat
      pure subroutine ten2mat_r4_DP(matrix, tensor, factors)
      use femtypes
      implicit none
      real(DP) :: matrix(6,6), tensor(3,3,3,3), factors(3)
      intent (in) :: tensor, factors
      intent (out) :: matrix
      end

      pure subroutine ten2mat_r4_DPC(matrix, tensor, factors)
      use femtypes
      implicit none
      real(DPC) :: factors(3)
      complex(DPC) :: matrix(6,6), tensor(3,3,3,3)
      intent (in) :: tensor, factors
      intent (out) :: matrix
      end
end interface


interface
      subroutine tensorgrid(px,py,xp,yp,maxlen,vx,vy,vxy,acolor)
      use femtypes
      implicit none
      integer (I4B) :: px, py, acolor
      real (DP) :: maxlen
      real (DP) :: xp(:,:), yp(:,:), vx(:,:), vy(:,:), vxy(:,:)
      intent (in) :: px, py, xp, yp, maxlen, vx, vy, acolor
      end subroutine tensorgrid
end interface


interface
      subroutine testac(anzzwg,anzknt,xpoint,ypoint,zki,lzrb,layrb)
      use femtypes
      implicit none
      integer (I4B) :: anzzwg, anzknt, zki(:,:), layrb(:,:), lzrb(:)
      real (DP) :: xpoint(:), ypoint(:)
      intent (in) :: anzzwg, anzknt, xpoint, ypoint, zki, lzrb, layrb
      end subroutine testac
end interface


interface
      subroutine text(accur,scalfk,anzgeb,xgeb,ygeb,txtgeb,layanz,      &
     &  lplayt,laytxt,layrb,txtlen)
      use femtypes
      implicit none
      integer (I4B) :: anzgeb, txtlen, layanz
      integer (I4B), pointer :: layrb(:,:), lplayt(:)
      real (DP) :: accur, scalfk
      real (DP), pointer :: xgeb(:), ygeb(:)
      character (len=txtlen), pointer :: txtgeb(:), laytxt(:)
      intent (in) :: accur, scalfk, txtlen
      intent (out) :: lplayt
      intent (inout) :: xgeb, ygeb, txtgeb
      intent (inout) :: anzgeb, layanz, laytxt, layrb
      end subroutine text
end interface


interface
      real (DP) function timestamp()
      use femtypes
      implicit none
      end function timestamp
end interface


interface
      subroutine total_loss(region,loss)
      use femtypes
      implicit none
      character (len=*), intent (in) :: region
      real (DP), intent (out) :: loss
      end subroutine total_loss
end interface


! TODO: replace old tensor routines by new ones
!interface
!      function trans_matrix(x1,x2)
!      use femtypes
!      implicit none
!      real (DP) :: x1(3), x2(3), trans_matrix(3,3)
!      intent (in) :: x1, x2
!      end function trans_matrix
!end interface


interface trans_matrix
      function trans_matrix(x1,x2)
      implicit none
      real (8) :: x1(3), x2(3), trans_matrix(3,3)
      intent (in) :: x1, x2
      end
end interface


interface transftensor
      pure subroutine transftensor_3Drank1(tensor1, transmat)
      use femtypes
      implicit none
      real (DP) :: tensor1(3)
      real (DP) :: transmat(3,3)
      intent (in) :: transmat
      intent (inout) :: tensor1
      end subroutine transftensor_3Drank1

      pure subroutine transftensor_3Drank2(tensor1, transmat)
      use femtypes
      implicit none
      real (DP) :: tensor1(3,3)
      real (DP) :: transmat(3,3)
      intent (in) :: transmat
      intent (inout) :: tensor1
      end subroutine transftensor_3Drank2
      
      pure subroutine transftensor_3Drank3(tensor1, transmat)
      use femtypes
      implicit none
      real (DP) :: tensor1(3,3,3)
      real (DP) :: transmat(3,3)
      intent (in) :: transmat
      intent (inout) :: tensor1
      end subroutine transftensor_3Drank3

      pure subroutine transftensor_3Drank4(tensor1, transmat)
      use femtypes
      implicit none
      real (DP) :: tensor1(3,3,3,3)
      real (DP) :: transmat(3,3)
      intent (in) :: transmat
      intent (inout) :: tensor1
      end subroutine transftensor_3Drank4
end interface


interface transformtensor
      pure subroutine transformtensor4_DP(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DP) :: a(3,3)
      real(DP) :: intensor(3,3,3,3), outtensor(3,3,3,3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
      end

      pure subroutine transformtensor3_DP(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DP) :: a(3,3)
      real(DP) :: intensor(3,3,3), outtensor(3,3,3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
      end

      pure subroutine transformtensor2_DP(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DP) :: a(3,3)
      real(DP) :: intensor(3,3), outtensor(3,3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
      end

      pure subroutine transformtensor1_DP(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DP) :: a(3,3)
      real(DP) :: intensor(3), outtensor(3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
      end

      pure subroutine transformtensor4_DPC(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DPC) :: a(3,3)
      complex(DPC) :: intensor(3,3,3,3), outtensor(3,3,3,3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
      end

      pure subroutine transformtensor3_DPC(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DPC) :: a(3,3)
      complex(DPC) :: intensor(3,3,3), outtensor(3,3,3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
      end

      pure subroutine transformtensor2_DPC(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DPC) :: a(3,3)
      complex(DPC) :: intensor(3,3), outtensor(3,3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
      end

      pure subroutine transformtensor1_DPC(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DPC) :: a(3,3)
      complex(DPC) :: intensor(3), outtensor(3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
      end
end interface transformtensor


interface
      subroutine transhom(a,tx,ty,tz)
      use femtypes
      implicit none
      real (DP) a(4,4), tx, ty, tz
      intent (in) tx, ty, tz
      intent (inout) a
      end subroutine transhom
end interface


interface
      subroutine transientsolver(acsr,acsrs,acsrd,acsrm,rhs,ia,ja,d2xdt2_at_tminusdt,stepcount,returncount, &
      &          trangtol,tranltollo,tranltolup,currenttime,timestep,timestepreject,epsgl,resgl)
      use femtypes
      implicit none
      complex (DPC), pointer :: acsr(:), acsrs(:), acsrd(:), acsrm(:), rhs(:)
      complex (DPC), pointer :: d2xdt2_at_tminusdt(:)
      integer (I4B), pointer :: ia(:), ja(:)
      integer (I4B) :: stepcount, returncount
      real (DP) :: trangtol, tranltollo, tranltolup, currenttime, timestep
      real (DP) :: epsgl, resgl
      logical :: timestepreject
      intent (in) :: trangtol
      intent (out) :: epsgl, resgl
      intent (inout) :: stepcount, returncount, currenttime, timestep, timestepreject
      intent (inout) ::  tranltollo, tranltolup
      end subroutine transientsolver
end interface


interface
      subroutine transsol(xnold,ynold,eold,enold,epold,egold,xold,nold)
      use femtypes
      implicit none
      integer (I4B) nold
      integer (I4B) :: eold(:,:), enold(:,:)
      integer (I4B) :: epold(:,:)
      real (DP) :: xnold(:), ynold(:)
      complex (DPC) :: xold(:)
      type (ARRPTRI) :: egold(:,:)
      intent (in) :: xnold, ynold, eold, enold, epold, egold, xold, nold
      end subroutine transsol
end interface


interface transten_r1
      pure function transten_r1_dp(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3), b(3), transten_r1_dp(3)
      intent (in) :: a, b
      end function transten_r1_dp
end interface


interface transten_r2
      pure function transten_r2_dpc(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3)
      complex (DPC) :: b(3,3), transten_r2_dpc(3,3)
      intent (in) :: a, b
      end function transten_r2_dpc

      pure function transten_r2_dp(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3), b(3,3), transten_r2_dp(3,3)
      intent (in) :: a, b
      end function transten_r2_dp
end interface


interface transten_r3
      pure function transten_r3_dpc(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3)
      complex (DPC) :: b(3,3,3), transten_r3_dpc(3,3,3)
      intent (in) :: a, b
      end function transten_r3_dpc

      pure function transten_r3_dp(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3), b(3,3,3), transten_r3_dp(3,3,3)
      intent (in) :: a, b
      end function transten_r3_dp
end interface


interface transten_r4
      pure function transten_r4_dpc(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3)
      complex (DPC) :: b(3,3,3,3), transten_r4_dpc(3,3,3,3)
      intent (in) :: a, b
      end function transten_r4_dpc

      pure function transten_r4_dp(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3), b(3,3,3,3), transten_r4_dp(3,3,3,3)
      intent (in) :: a, b
      end function transten_r4_dp
end interface


interface
      subroutine tria(xk,yk,plist,plp,zahl,prlist,pzr,prz,ele,ez,       &
     &  ezalt,kdim,geb,nb,bz,erfolg)
      use femtypes
      implicit none
      integer (I4B) plist(:),plp(:),zahl,prlist(:,:),pzr,prz
      integer (I4B) ele(:,:),ez,ezalt,kdim,geb(:),nb(:,:),bz
      real (DP) xk(:),yk(:)
      logical erfolg
      intent (in) :: xk, yk, zahl, pzr, ezalt, kdim, bz
      intent (out) :: erfolg
      intent (inout) :: plist, plp, prlist, prz, ele, ez, geb, nb
      end subroutine tria
end interface


interface
      subroutine tripc(zr,zc,nc,x,y,e,en,ie)
      use femtypes
      implicit none
      real (DP) x(:), y(:), zr(:), zc(:)
      integer (I4B) e(:,:), en(:,:), ie, nc
      intent (in):: zr, zc, nc, x, y, e, en, ie
      end subroutine tripc
end interface


interface
      subroutine triplt(e,en,x,y,ie,r,g,b)
      use femtypes
      implicit none
      real (DP) x(:), y(:)
      real (SP) r, g, b
      integer (I4B) e(:,:), en(:,:), ie
      intent (in)::e, en, x, y, ie, r, g, b
      end subroutine triplt
end interface


interface umfcsolver
      subroutine UMFCSOLVER_Z(a,b,x,n,ia,ja)
      use femtypes
      implicit none
      complex(DPC), pointer :: a(:), b(:)
      complex(DPC), pointer :: x(:)
      integer(I4B), pointer :: ia(:),ja(:)
      integer(I4B) n
      intent(in) :: n
      intent(out) :: x
      intent(inout) :: a, b, ia, ja
      end subroutine UMFCSOLVER_Z

      subroutine UMFCSOLVER_D(a,b,x,n,ia,ja)
      use femtypes
      implicit none
      real(DP), pointer :: a(:), b(:)
      real(DP), pointer :: x(:)
      integer(I4B), pointer :: ia(:),ja(:)
      integer(I4B) n
      intent(in) :: n
      intent(out) :: x
      intent(inout) :: a, b, ia, ja
      end subroutine UMFCSOLVER_D
end interface


interface umfsolver
      subroutine umfsolver_z(a,b,x,n,ir,ic)
      use femtypes
      implicit none
      complex (DPC), pointer :: a(:)
      complex (DPC) b(:),x(:)
      integer (I4B) n
      integer (I4B), pointer :: ir(:),ic(:)
      intent (in) :: n
      intent (inout) :: b, x, ir, ic, a
      end subroutine umfsolver_z

      subroutine umfsolver_d(a,b,x,n,ir,ic)
      use femtypes
      implicit none
      real (DP), pointer :: a(:)
      real (DP) :: b(:),x(:)
      integer (I4B) :: n
      integer (I4B), pointer :: ir(:),ic(:)
      intent (in) :: n
      intent (inout) :: b, x, ir, ic, a
      end subroutine umfsolver_d

      subroutine umfsolver_s(a,b,x,n,ir,ic)
      use femtypes
      implicit none
      real (SP), pointer :: a(:)
      real (SP) :: b(:),x(:)
      integer (I4B) :: n
      integer (I4B), pointer :: ir(:),ic(:)
      intent (in) :: n
      intent (inout) :: b, x, ir, ic, a
      end subroutine umfsolver_s
end interface umfsolver


interface umx
      subroutine umx_c(lower,ia,ja,v1,e,n,om1)
      use femtypes
      implicit none
      complex (SPC) lower(:),v1(:),e(:)
      real (SP) om1
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, ia, ja, v1, n, om1
      intent (out) :: e
      end subroutine umx_c

      subroutine umx_z(lower,ia,ja,v1,e,n,om1)
      use femtypes
      implicit none
      complex (DPC) lower(:),v1(:),e(:)
      real (DP) om1
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, ia, ja, v1, n, om1
      intent (out) :: e
      end subroutine umx_z

      subroutine umx_d(lower,ia,ja,v1,e,n,om1)
      use femtypes
      implicit none
      integer (I4B) ia(:),ja(:),n
      real (DP) lower(:),v1(:),e(:),om1
      intent (in) :: lower, ia, ja, v1, n, om1
      intent (out) :: e
      end subroutine umx_d

      subroutine umx_s(lower,ia,ja,v1,e,n,om1)
      use femtypes
      implicit none
      integer (I4B) ia(:),ja(:),n
      real (SP) lower(:),v1(:),e(:),om1
      intent (in) :: lower, ia, ja, v1, n, om1
      intent (out) :: e
      end subroutine umx_s
end interface


interface
      subroutine unread
      end subroutine unread
end interface


interface unscal
      subroutine unscal_c(lower,upper,diag,b,x,ia,ja,n,resnrm,bscal,h)
      use femtypes
      implicit none
      complex (SPC) lower(:),upper(:),diag(:),x(:),b(:),h(:)
      real (SP) resnrm,bscal
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: lower, upper, diag, ia, ja, n, bscal
      intent (out) :: resnrm, h
      intent (inout) ::  b, x
      end subroutine unscal_c

      subroutine unscal_z(lower,upper,diag,b,x,ia,ja,n,resnrm,bscal,h)
      use femtypes
      implicit none
      complex (DPC) lower(:),upper(:),diag(:),x(:),b(:),h(:)
      real (DP) resnrm,bscal
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: lower, upper, diag, ia, ja, n, bscal
      intent (out) :: resnrm, h
      intent (inout) ::  b, x
      end subroutine unscal_z

      subroutine unscal_d(lower,upper,diag,b,x,ia,ja,n,resnrm,bscal,h)
      use femtypes
      implicit none
      real (DP) lower(:),upper(:),diag(:),x(:),b(:), bscal
      real (DP) h(:)
      real (DP) resnrm
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: lower, upper, diag, ia, ja, n, bscal
      intent (out) :: resnrm, h
      intent (inout) ::  b, x
      end subroutine unscal_d

      subroutine unscal_s(lower,upper,diag,b,x,ia,ja,n,resnrm,bscal,h)
      use femtypes
      implicit none
      real (SP) lower(:),upper(:),diag(:),x(:),b(:), bscal
      real (SP) h(:)
      real (SP) resnrm
      integer (I4B) n,ia(:),ja(:)
      intent (in) :: lower, upper, diag, ia, ja, n, bscal
      intent (out) :: resnrm, h
      intent (inout) ::  b, x
      end subroutine unscal_s
end interface


interface
      subroutine unsect
      use femtypes
      implicit none
      end subroutine unsect
end interface


interface
      subroutine userbc(branch,bcnum,xs,ys,nature,pval,qval,elem)
      use femtypes
      implicit none
      integer (I4B) branch, bcnum, nature
      integer (I4B) , optional :: elem
      real (DP) xs, ys
      complex (DPC) pval
      complex (DPC), optional :: qval(:)
      intent (in) :: branch, bcnum, xs, ys, nature, elem
      intent (out) :: pval, qval
      end subroutine userbc
end interface


interface usermaterials
      pure subroutine usermaterials(matname,found,xys,names,values,numnames,elem)
      use femtypes
      implicit none
      real (DP), optional :: xys(:), values(:)
      integer (I4B),optional :: elem
      integer (I4B), optional :: numnames
      character (len=*) :: matname
      character (len=*), optional :: names(:)
      logical :: found
      intent (in) :: elem, matname, xys
      intent (out) :: names, values, found, numnames
      end subroutine usermaterials
end interface


interface
      subroutine vertex(x2,y2,z2,ausbuc,scalfk,layanz,                  &
     &  laytxt,layrb,layer,txtlen)
      use femtypes
      implicit none
      integer (I4B) :: layanz, layer, txtlen
      integer (I4B), pointer :: layrb(:,:)
      real (DP) :: x2, y2, z2, ausbuc, scalfk
      character (len=txtlen), pointer :: laytxt(:)
      intent (in) :: scalfk, txtlen
      intent (out) :: x2, y2, z2, ausbuc, layer
      intent (inout) :: layanz, laytxt, layrb
      end subroutine vertex
end interface


interface
      subroutine wfein(ellist,zahle,zahlk,maxelemt,ende)
      use femtypes
      implicit none
      integer (I4B) ellist(:), zahle, zahlk, maxelemt
      logical ende
      intent (in) :: maxelemt
      intent (out) :: ellist, ende
      intent (inout) :: zahle, zahlk
      end subroutine wfein
end interface


interface
      pure real (DP) function wink(x,y,x1,y1,x2,y2)
      use femtypes
      implicit none
      real (DP) x,y,x1,y1,x2,y2
      intent (in) :: x, y, x1, y1, x2, y2
      end function wink
end interface


interface
      pure real (DP) function wink1(x,y,x1,y1,x2,y2)
      use femtypes
      implicit none
      real (DP) x,y,x1,y1,x2,y2
      intent (in) :: x, y, x1, y1, x2, y2
      end function wink1
end interface


interface
      real (DP) function winkl(la,mi,ne,xn,yn)
      use femtypes
      implicit none
      integer (I4B) la, mi, ne
      real (DP) xn(:), yn(:)
      intent (in) :: la, mi, ne, xn, yn
      end function winkl
end interface


interface
      subroutine wnetin(gbz,gzz,gkz,bzi,bzip,bzil,zki,kzrb,zpz,zpp,     &
     &  xbk,ybk,matname,filenm,lzrb,layrb,alrb,btrb,laytxt,layanz,nnat)
      use femtypes
      implicit none
      integer (I4B) :: gbz, gzz, gkz, nnat, layanz
      integer (I4B) :: bzi(:), bzip(:), bzil(:), zki(:,:), kzrb(:,:)
      integer (I4B) :: zpz(:), lzrb(:), layrb(:,:)
      real (DP) :: xbk(:), ybk(:), zpp(:)
      character (len=*) :: filenm, matname(:), laytxt(:)
      complex (DPC) :: alrb(:,:), btrb(:,:)
      intent (in) :: gbz, gzz, gkz, bzi, bzip, bzil, zki, kzrb, zpz, zpp
      intent (in) :: xbk, ybk, matname, filenm, lzrb, layrb, alrb, btrb
      intent (in) :: laytxt, layanz, nnat
      end subroutine wnetin
end interface


interface
      subroutine wpathsetting
      use femtypes
      implicit none
      end subroutine wpathsetting
end interface


interface
      subroutine writsetting
      use femtypes
      implicit none
      end subroutine writsetting
end interface


interface
      pure subroutine xy2lam(xt,yt,elem,lambda,xn,yn,e)
      use femtypes
      implicit none
      integer (I4B) elem, e(:,:)
      real (DP) xt, yt, lambda(3), xn(:), yn(:)
      intent (in) :: xt, yt, xn, yn, elem, e
      intent (out) :: lambda
      end subroutine xy2lam
end interface


interface
      subroutine xyzout(ie,e,x,y,p,z,is,winkel)
      use femtypes
      implicit none
      integer (I4B) p,ie,is,e(:,:)
      real (DP) winkel,x(:),y(:)
      complex (DPC) z(:)
      intent (in) :: ie ,e ,x ,y, p ,z, is, winkel
      end subroutine xyzout
end interface


interface
      subroutine zanfpp(plottyp,xmin,xmax,ymin,ymax,h,bildnr,information)
      use femtypes
      implicit none
      real (DP) xmin,xmax,ymin,ymax,h
      integer (I4B) bildnr
      character (len=*) plottyp
      character (len=*),optional :: information
      intent (in):: plottyp, information
      intent (inout):: xmin,xmax,ymin,ymax
      intent (out):: bildnr,h
      end subroutine zanfpp
end interface


interface
      subroutine zeit(text)
      use femtypes
      implicit none
      character (len=*) text
      intent (in) :: text
      end subroutine zeit
end interface


interface zqmrsolver
      subroutine zqmrsolver_DPC(a,b,xzr,n,eps,ia,ja,epsgl,resgl,symm)
      use femtypes
      implicit none
      complex (DPC) a(:),b(:),xzr(:)
      integer (I4B) n,ia(:),ja(:)
      real (DP) eps,epsgl,resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: a, b, xzr
      end subroutine zqmrsolver_DPC

      subroutine zqmrsolver_DP(a,b,xzr,n,eps,ia,ja,epsgl,resgl,symm)
      use femtypes
      implicit none
      real (DP) a(:),b(:),xzr(:)
      integer (I4B) n,ia(:),ja(:)
      real (DP) eps,epsgl,resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: a, b, xzr
      end subroutine zqmrsolver_DP
end interface

end module feminterface
