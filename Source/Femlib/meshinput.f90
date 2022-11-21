      subroutine meshinput(ok,jd,reell,gilt,istlin,konfom,zyl,epsgl,resgl,resnl,error)
  
      use feminterface, only: readnetin, extend, basin
      use femtypes
      use globalvariables, only: gbz,gzz,gkz,zki,kzrb,zrb,xbk,ybk
      use globalvariables, only: matzif,alrb,btrb,nnat,xmin,xmax,ymin,ymax
      use globalvariables, only: n,e,p,xn,yn,kzi,geb,en
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
!    $Revision: 1.13 $
!    $Date: 2014/07/15 13:04:51 $
!    $Author: m_kasper $
!
! meshinput   Reads mesh information from basis

      real (DP) epsgl, resgl, resnl, error, jd
      logical ok, reell, gilt, istlin, konfom, zyl

      intent (out) :: ok, reell, gilt, istlin, konfom, zyl, jd
      intent (out) :: epsgl, resgl, resnl, error

!===================================================================================
! ELEMENTS:
!
!     e           element information, nodes of the elements
!     en          neighbors of the elements
!     geb         Gebiets (area) number for each element
!     n           number of the elements
!---------------------------------------------------------------------------------
! NODES:
! 
!     xn,yn       nodes coordinates
!     kzi         Knoten Zweig Information (node branch information)
!                      /   = 0  : inner nodes
!                 kzi  -   < 0  : node is identical to keypoint (keypoint index)
!                      \   > 0  : node is on branch. (kzi)
!     p           number of the nodes
!----------------------------------------------------------------------------------
! SOLUTION:
!
!     ep          Polynomial degree of elements ( TO DO)
!     eg          Global degrees of freedom for the elements (vector of vectors)
!                 eg(i)%d(j) is a mapping between the local index j of element i and
!                 the global dof index
!     x           Solution vector ( TO DO in Solver)
!     ndof        Total number of degrees of freedom
!---------------------------------------------------------------------------------
! BRANCHES:
!
!     zki         Zweige-nodes-Informationen (Liste the nodes zu jedem Zweig)
!                 zki(1,i)  Anfangspunkt the Zweiges
!                 zki(2,i)  Endpunkt the Zweiges
!                 zki(3,i)  Dritter Punkt the Zweiges
!                           /  = 0  :  Gerade
!                 zki(3,i)  -  > 0  :  Punkt auf dem Kreisbogen
!                           \  < 0  :  Mittelpunkt the Kreisbogens
!     zrb         Type of BC for the branch:
!                     0 -  99: Dirichlet
!                   200 - 299: Neumannn
!                   300 - 399: Kontur
!                   400 - 499: innen
!     alrb, btrb  (complex) parameter of the Dirichlet of general BC
!     gzz         Gebiets-Zweige-Zahl (number of branches)
!--------------------------------------------------------------------------------
! KEYPOINTS:
!
!     xbk,ybk     coordinates of keypoints
!     kzrb        boundary condition of the keypoints, index of the branch whose BC is
!                 valid for the keypoint
!     gkz         Gebiets-Knoten-Zahl (number of nodes)
!-----------------------------------------------------------------------------------
! REGIONS:
!
!     matzif      material index with which the region is filled
!     matname     (to implement)
!     gbz         Gebiets-Bereiche-Zahl (number of regions)
!     bzi         Bereichs-Zweige-Informationen (list of branches which form the regions)
!     bzip        pointer to the start of branch of a region in bzi (for compact storage)
!     bzil        pointer onto bzip including all regions inside a region (multiply
!                 connected regions )
!-------------------------------------------------------------------------------------
! ETC:
!
!     ok          =.false. if an error is found
!     omega       radial frequency
!     maxgbz      maximum sizes for gbz
!     maxgzz      maximum sizes for gzz
!     maxgkz      maximum sizes for gkz
!

!     reell       =.true. if the solution vector is real (omega=0)
!     gilt        =.true. if the file "solution" contains a valid solution vector
!     istlin      =.true. if contains no nonlinear material properties 
!     konform     =.true. if conforming elements are used
!     zyl         =.true. if it handles a  cylindrical symmetric problem
!     epsgl       solution precision of the linear equation systems 
!                 (precision of the solution vector)
!     resgl       residual of the equation system
!     resnl       residual of the nonlinear equation (nonlinear materials)
!     error       bezogener average error of the solution 
!                 (from the aposteriori error abschaetzung)
!------------------------------------------------------------------------------------
! AUXILIARY:
!
!     xmin, xmax  aeusserste Ausdehnungen the Gebiets
!     ymin, ymax  
!====================================================================================

!  lokal variables
      real (DP) omegaf
      logical ok2
!
!  read the geometry data: key-points, branches, regions
!
      call readnetin(gbz,gzz,gkz,zki,kzrb,zrb,xbk,ybk,matzif,ok2,alrb,btrb,nnat)
!  compute the extends (x-, y-bounds of geometry)
!
      call extend(gzz,xbk,ybk,zki,xmin,xmax,ymin,ymax)
!  read mesh data: elements nodes
!
      call basin(ok,jd,n,p,e,xn,yn,kzi,geb,en)
      write (*,'(a,i7)') '  mesh nodes:',p
      write (*,'(a,i7)') '  elements:  ',n
!
!  check data
      if ( (.not. ok2) .or. (.not. ok) .or. (n.le.0) .or. (p.le.0) ) then
        print*,'**** error in the mesh data'
        ok=.false.
        stop
      end if 
!
      reell=.false.
      omegaf=0._DP
      istlin=.false.
      gilt=.false.
      konfom=.false.
      zyl=.false.
      epsgl=1._DP
      resgl=1._DP
      resnl=1._DP
      error=huge(1._DP)
      return
      end subroutine meshinput
      
