      subroutine nfarbep(palette,colorfill,ncolor,zcolmin,zcolmax,      &
     &  hlin,nlin,zlinmin,zlinmax,r,g,b,phi,fieldtype)
      use feminterface, only : drtri, fieldquantity, hline
      use globalvariables
      use femtypes
      implicit none
      real (SP) :: r, g, b
      real (DP) :: zcolmin, zcolmax, zlinmin, zlinmax, phi
      integer (I4B) :: ncolor, nlin, palette(:)
      logical :: colorfill, hlin
      character (len=*) :: fieldtype
      intent(in):: palette, colorfill, ncolor, zcolmin, zcolmax, hlin, nlin
      intent(in):: zlinmin, zlinmax, r, g, b, phi, fieldtype
!
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
!    $Revision: 1.35 $
!    $Date: 2014/07/01 14:43:50 $
!    $Author: m_kasper $
!
!------------------------------------------------------------------------------
!
!  Color fill and equipotentials on a triangle mesh   2D Plot
!  Coloring according to the potential values  z  (solution vector)
!  Subdivision of triangles into subtriangles and linear interpolation of the color
!
!       z           potential
!       xn          vertex x-coordinates of triangles
!       yn          vertex x-coordinates of triangles
!       e           vertex numbers of triangles
!       n           number of triangles 
!       palette     color palette in increasing order
!       colorfill   =.true, for color fill plot
!       ncolor      number of colors
!       zcolmin     potential of the lowest color
!       zcolmax     potential of the largest color
!       hlin        =.true. if equipotential should be plotted
!       nlin        number of equipotentials
!       zlinmin     potential of the lowest equipotential
!       zlinmax     potential of the largest equipotential
!       r, g, b     red green and blue value for the potential lines
!       phi         time point (t*omega) in degree
!       fieldtype   defines which kind of field quantity should be plotted
!
      integer (I4B) i, oldindex, elelev
      integer (I4B) is, lev, maxlev, k, k1, k2, arrsize
      integer (I4B), allocatable:: es(:,:) 
      real (SP) ocr, ocg, ocb, olcr, olcg, olcb
      real (DP), allocatable :: zcolor(:), zhlin(:)
      real (DP), allocatable:: xs(:), ys(:), zs(:)
      real (DP) ux, uy, vx, vy, lam(3)
      complex (DPC) :: fieldval
      logical proj, ok1
      integer (I4B), parameter :: maxplevel=20
      integer (I4B) ,parameter :: plev(maxplevel)=(/1,5,9,12,15,17,19,21,22,23,24,25,26,27,28,29,30,31,32,33/)
      character (len=10) :: unit
      character (len=50) :: descriptor
!
      external :: grglun, pgbeg, pgbox, pgcurs, pgdraw, pgend, pglabel, pgline, pgsch
      external :: pgqci, pgscf, pgsci, pgscr, pgsfs, pgslw, pgqch, pgqcr, pgmove
      external :: pgopen, pgpap, pgqcs, pgqwin, pgqvp, pgsvp, pgswin, pgwnad
      external :: pgrect, pgptxt

      if (.not. (colorfill.or.hlin) ) return
      proj= .false.
      call pgqci(oldindex)
      if (colorfill) then 
!  save color status, use color index 20 for filling
        call pgsci(20)
        call pgqcr(20, ocr, ocg, ocb)
!  set filling style
        call pgsfs(1)
!  assignment of colors to the potential values
        allocate(zcolor(ncolor))
        do i=1,ncolor
          zcolor(i)=(zcolmax-zcolmin)/real(ncolor,DP)*real(i,DP)+zcolmin
        end do
      end if
      if (hlin) then 
!  save color status, use color index 21 for equipotentials
        call pgsci(21)
        call pgqcr(21, olcr, olcg, olcb)
        call pgscr(21,  r, g, b)
!  assignment of the potential values
        allocate(zhlin(nlin))
        do i=1,nlin
          zhlin(i)=(zlinmax-zlinmin)/real(nlin-1,DP)*real(i-1,DP)+zlinmin
        end do
      end if
!  compute subelements
      maxlev=plev(maxplevel)
      allocate(es(3,maxlev**2))
      do lev=1,maxlev
        do is=(lev-1)**2+1,lev**2
          es(1,is)=int(real(is+lev)/2.)
          es(2,is)=int(real(is+1+lev)/2.)+lev
          es(3,is)=int(real(is+lev)/2.)+1+mod(is+lev+1,2)*lev
        end do
      end do
      arrsize=(maxlev+1)*(maxlev+2)/2
      allocate(xs(arrsize),ys(arrsize),zs(arrsize))
!  loop over all elements
      do i=1,n
!  set the level of subdivision in dependence of the polynomial degree
        elelev=plev(min(maxval(ep(i,:)),maxplevel))
!  compute the subelement coordinates
        ux=( xn(e(2,i))-xn(e(1,i)) ) / real(elelev,DP)
        uy=( yn(e(2,i))-yn(e(1,i)) ) / real(elelev,DP)
        vx=( xn(e(3,i))-xn(e(2,i)) ) / real(elelev,DP)
        vy=( yn(e(3,i))-yn(e(2,i)) ) / real(elelev,DP)
!
!$omp parallel do default(none) &
!$omp shared( elelev, xs, ys, zs, ux, uy, vx, vy, e, i) &
!$omp shared( fieldtype, phi, physics, xn, yn ) &
!$omp private( lam, fieldval, descriptor, unit, ok1, k2, k1)
        do k=1,(elelev+1)*(elelev+2)/2
          k1=int(sqrt(real(2*k))-0.5)
          k2= k - k1*(k1+1)/2 -1
          xs(k)=k1*ux+k2*vx+xn(e(1,i))
          ys(k)=k1*uy+k2*vy+yn(e(1,i))
          lam(1)=1._DP-real(k1,DP)/real(elelev,DP)
          lam(3)=real(k2,DP)/real(elelev,DP)
          lam(2)=1._DP-lam(1)-lam(3)
!
          call fieldquantity(i,fieldtype,lam,phi,fieldval,descriptor,unit,ok1)
          if (.not.ok1) then
            print *,"***** fieldtype: ",trim(fieldtype)," unknown for ",trim(physics)
            stop
          end if
          zs(k) = real(fieldval,DP)
        end do
!$omp end parallel do 
!
!  loop over all subelements: color fill
        if (colorfill) then
          call pgsci(20)
          do is=1,elelev**2          
            call drtri(is,es,xs,ys,zs,zcolor,ncolor,palette,proj)
          end do
        end if
!  loop over all subelements: equipotentials 
        if (hlin) then
          call pgsci(21)
          do is=1,elelev**2          
            call hline(is,es,xs,ys,zs,zhlin,nlin,proj)
          end do
        end if
      end do
      deallocate(xs, ys, zs, es)
!  reset color status
      if (colorfill) then 
        deallocate(zcolor)
        call pgscr(20, ocr, ocg, ocb)
      end if
      if (hlin) then 
        deallocate(zhlin)
        call pgscr(21, olcr, olcg, olcb)
      end if
      call pgsci(oldindex)
      return
      end subroutine nfarbep