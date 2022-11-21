      program triang
!
!  program to generate a triangulation
!
      use feminterface, only: initialize, getsetting, zeit, lese,       &
     &  sorttr, wnetin, minmesh, ausgab, delaunay_refinement, glaett,   &
     &  zanfpp, struc, triplt, farbpalette, farbe, angles, extend, nachb, dichk
      use globalvariables
      use femtypes
      implicit none
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
!    $Revision: 1.24 $
!    $Date: 2010/09/06 12:41:35 $
!    $Author: m_kasper $
!
      integer (I4B) :: bildnr, layanz, fehlkn, i, ncolor
      integer  (I4B), allocatable:: palette(:), liste5(:)
      integer (I4B), pointer :: bzi(:), bzip(:), bzil(:), krb(:,:)
      integer (I4B), pointer :: zpz(:), lzrb(:), layrb(:,:)
      real (DP) :: rand, xpl, xpu, ypl, ypu, w1, w2, w3
      real (DP) :: h, zcolmin, zcolmax
      real (SP) :: r = 0.00, g = 0.00, b = 0.00
      real (DP), allocatable :: wmin(:)
      real (DP), pointer :: zpp(:), meshsize(:)
      complex (DPC), pointer :: lalrb(:,:), lbtrb(:,:)
      logical :: ok, konfom, zyl, nregen
      character (len=25) :: meshsmooth, plotdevice
      character (len=25), pointer ::  matname(:), laytxt(:)
      character (len=200) :: path
!
!  explanation of some variables:
!
!  Regions are stored in a compact format using three arrays: bzil, bzip, bzi
!    gbz        number of regions
!    bzi        is a continous list of branches, including the inner branches. 
!               The sequence is the same as in the input file (netin.acd)
!    bzip       is a pointer (index) to the first branch in the list bzi of a 
!               region. For multiply connected regions inner boundary branches 
!               precede the outer boundary.  
!    bzil       for each region is a pointer (index) to the first 
!               (possibly inner) region in the in the array bzip
!    matname    materialnames of the regions
!  Branches
!    gzz        number of branches
!    zki        key-points (geometry nodes) which determine the branches
!               zki(1,i)  start key-point of the branch
!               zki(2,i)  end key-point of the branch
!               zki(3,i)  third key-point of the branch, with
!                         /  = 0  :  straight line 
!               zki(3,i)  -  > 0  :  node located on an arc              
!                         \  < 0  :  midpoint of an arc
!    zrb        type of the boundary condition of the branch, 
!               which can take the following values:
!                 0 -  99: Dirichlet BC, where for 1- 99: the BC has to be set in the file USERBC 
!               200 - 299: Neumann of general BC
!               300 - 399: visible contour, no BC
!               400 - 499: non-visible contour, no BC
!    zpp        control parameter for the distribution of nodes
!                    /  = 0  :  equidistant distribution
!               zpp  -  > 0  :  nodes are denser at the startpoint
!                    \  < 0  :  nodes are denser at the endpoint
!    zpz        number of nodes on the branch (including start- and endpoint)
!  Key-points
!    gkz        number of key-points
!    xbk,ybk    coordinates of key-points
!    krb        boundary condition of key-points
!    kzrb       branch number which determines the boundary condition of a key-point
!  Nodes
!    kzi        node branch information
!               kzi=0: the node is an inner node
!               kzi<0: the node is identical to the keypoint with number (-kzi)
!               kzi>0: the node is on the branch with the number (kzi)
!    xn, yn     coordinates of the nodes
!  Auxillary
!    xmin,xmax  extends of the problem,
!    ymin,ymax  outermost bounds of the geometry in x- and y-direction
!
      call initialize()
!
      call getsetting('PROJECTPATH',path)
      ok=.true.
!  read input data for the triangulation
      call lese(bzi,bzip,bzil,zpz,zpp,matname,krb,meshsize,             &
     &                lzrb,layanz,layrb,laytxt,lalrb,lbtrb)
!  correct if necessary direction of branches
      call sorttr(gbz,gzz,bzi,bzip,bzil,zki,zrb,ok)
!  write structure of the problem to the file netin.dat
      call wnetin(gbz,gzz,gkz,bzi,bzip,bzil,zki,kzrb,zpz,zpp,xbk,ybk,   &
     &            matname,path(1:len_trim(path))//'netin.dat',          &
     &            lzrb,layrb,lalrb,lbtrb,laytxt,layanz,nnat)
!
      call extend(gzz,xbk,ybk,zki,xmin,xmax,ymin,ymax)
!  margin in percent of the total size
      rand=.02_DP
      xpl=xmin-(xmax-xmin)*rand/2._DP
      xpu=xmax+(xmax-xmin)*rand/2._DP
      ypl=ymin-(ymax-ymin)*rand/2._DP
      ypu=ymax+(ymax-ymin)*rand/2._DP
      print * ,whatsystem
      select case (whatsystem)
      case ('WINDOWS')
        plotdevice='/WZ'
      case ('LINUX')
        plotdevice='/XWINDOW'
      case default
        plotdevice='/NULL'
      end select
      call zanfpp(plotdevice,xpl,xpu,ypl,ypu,h,bildnr)
!  draw the structure
      call struc(gzz,zki,zrb,xbk,ybk,r,g,b)
!
!  draw the structure with bridges
      call zanfpp(plotdevice,xpl,xpu,ypl,ypu,h,bildnr)
!  draw the structure
      call struc(gzz,zki,zrb,xbk,ybk,r,g,b)
!  generate the mesh
      call minmesh(bzi,bzip,bzil,zpz,zpp,meshsize)
!  draw initial minimal mesh
      call zanfpp(plotdevice,xpl,xpu,ypl,ypu,h,bildnr)
!
      allocate (en(3,size(geb)))     
      call nachb(e,en,p,n)
      call triplt(e,en,xn,yn,n,0.80,0.80,0.80)
      konfom=.false.
      zyl=.false.
!  mesh refinement
      call delaunay_refinement(meshsize)
!  mesh smoothing
      call getsetting('MESHSMOOTHING',meshsmooth)
      if (meshsmooth .eq. 'YES') then
        allocate (liste5(p), x(p))
        nregen=.false.
        fehlkn=0
        liste5(:)=0
        x(:)=(0.,0.)
        call glaett(n,e,xn,yn,en,geb,kzi,p,nregen,liste5,               &
     &    fehlkn,zrb,zki,xbk,ybk,20)
        deallocate(liste5,x)
      end if
!  sort, output mesh and startvector for solution
      call ausgab(e,xn,yn,en,n,p,kzi,zrb,krb,geb,konfom,zyl)
!  draw the mesh
      call zanfpp(plotdevice,xpl,xpu,ypl,ypu,h,bildnr)
      call triplt(e,en,xn,yn,n,0.80,0.80,0.80)
      call struc(gzz,zki,zrb,xbk,ybk,r,g,b)
!
!  color plot of mesh quality according to the minimum angle in the element
      call zanfpp(plotdevice,xpl,xpu,ypl,ypu,h,bildnr)
      allocate(wmin(n))
      do i=1,n
        call angles(i,e,xn,yn,w1,w2,w3)
        wmin(i)=min(w1,w2,w3)
      end do 
      ncolor=10
      allocate (palette(ncolor))
      call farbpalette(palette,ncolor,'CS')
!  blue is cool i.e. good mesh quality, red is hot small angle here
      zcolmin=3.1415926/3.
      zcolmax=0.
      call farbe(wmin,xn,yn,e,n,palette,ncolor,zcolmin,zcolmax)
      deallocate (palette, wmin)
      if (n < 20.e3) then
        call triplt(e,en,xn,yn,n,0.80,0.80,0.80)
      end if
      call struc(gzz,zki,zrb,xbk,ybk,r,g,b)
!
      call pgend
      stop
      end program triang
