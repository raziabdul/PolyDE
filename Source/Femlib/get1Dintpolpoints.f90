      pure subroutine get1Dintpolpoints(norder, xtab, errcode)
      use femtypes
      implicit none
      integer (I4B) norder, errcode
      real (DP) :: xtab(:)
      intent (in) :: norder
      intent (out) :: xtab, errcode
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
!    $Revision: 1.9 $
!    $Date: 2015/11/03 15:29:05 $
!    $Author: m_kasper $
!
!  This subroutine returns 1D Gauss-Lobatto integration points on the interval [-1,1]
!  used as INTERPOLATION POINTS.
!
!  Input:
!            norder   order of the rule (number of points and weights) in the range 2 to 22
!  Output:
!            xtab     vector of coordinates of the interpolation points 
!            errcode  =1001 if an integration routine of this order is not available
!  
!  Modified from routine lobatto_set of John Burkardt
!  http://www.math.iastate.edu/burkardt/f_src/quadrule/quadrule.html
!
!  References:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!
      select case (norder)
      case (1)
        xtab(1) =    0.0e+00_DP
      case (2)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =    1.0e+00_DP
      case (3)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =    0.0e+00_DP
        xtab(3) =    1.0e+00_DP
      case (4)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.447213595499957939281834733746e+00_DP
        xtab(3) =    0.447213595499957939281834733746e+00_DP
        xtab(4) =    1.0e+00_DP
      case (5)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.654653670707977143798292456247e+00_DP
        xtab(3) =    0.0e+00_DP
        xtab(4) =    0.654653670707977143798292456247e+00_DP
        xtab(5) =    1.0e+00_DP
      case (6)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.765055323929464692851002973959e+00_DP
        xtab(3) =  - 0.285231516480645096314150994041e+00_DP
        xtab(4) =    0.285231516480645096314150994041e+00_DP
        xtab(5) =    0.765055323929464692851002973959e+00_DP
        xtab(6) =    1.0e+00_DP
      case (7)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.830223896278566929872032213967e+00_DP
        xtab(3) =  - 0.468848793470714213803771881909e+00_DP
        xtab(4) =    0.0e+00_DP
        xtab(5) =    0.468848793470714213803771881909e+00_DP
        xtab(6) =    0.830223896278566929872032213967e+00_DP
        xtab(7) =    1.0e+00_DP
      case (8)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.871740148509606615337445761221e+00_DP
        xtab(3) =  - 0.591700181433142302144510731398e+00_DP
        xtab(4) =  - 0.209299217902478868768657260345e+00_DP
        xtab(5) =    0.209299217902478868768657260345e+00_DP
        xtab(6) =    0.591700181433142302144510731398e+00_DP
        xtab(7) =    0.871740148509606615337445761221e+00_DP
        xtab(8) =    1.0e+00_DP
      case (9)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.899757995411460157312345244418e+00_DP
        xtab(3) =  - 0.677186279510737753445885427091e+00_DP
        xtab(4) =  - 0.363117463826178158710752068709e+00_DP
        xtab(5) =    0.0e+00_DP
        xtab(6) =    0.363117463826178158710752068709e+00_DP
        xtab(7) =    0.677186279510737753445885427091e+00_DP
        xtab(8) =    0.899757995411460157312345244418e+00_DP
        xtab(9) =    1.0e+00_DP
      case (10)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.919533908166458813828932660822e+00_DP
        xtab(3) =  - 0.738773865105505075003106174860e+00_DP
        xtab(4) =  - 0.477924949810444495661175092731e+00_DP
        xtab(5) =  - 0.165278957666387024626219765958e+00_DP
        xtab(6) =    0.165278957666387024626219765958e+00_DP
        xtab(7) =    0.477924949810444495661175092731e+00_DP
        xtab(8) =    0.738773865105505075003106174860e+00_DP
        xtab(9) =    0.919533908166458813828932660822e+00_DP
        xtab(10) =   1.0e+00_DP
      case (11)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.934001430408059134332274136099e+00_DP
        xtab(3) =  - 0.784483473663144418622417816108e+00_DP
        xtab(4) =  - 0.565235326996205006470963969478e+00_DP
        xtab(5) =  - 0.295758135586939391431911515559e+00_DP
        xtab(6) =    0.0e+00_DP
        xtab(7) =    0.295758135586939391431911515559e+00_DP
        xtab(8) =    0.565235326996205006470963969478e+00_DP
        xtab(9) =    0.784483473663144418622417816108e+00_DP
        xtab(10) =   0.934001430408059134332274136099e+00_DP
        xtab(11) =   1.0e+00_DP
      case (12)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.944899272222882223407580138303e+00_DP
        xtab(3) =  - 0.819279321644006678348641581717e+00_DP
        xtab(4) =  - 0.632876153031869677662404854444e+00_DP
        xtab(5) =  - 0.399530940965348932264349791567e+00_DP
        xtab(6) =  - 0.136552932854927554864061855740e+00_DP
        xtab(7) =    0.136552932854927554864061855740e+00_DP
        xtab(8) =    0.399530940965348932264349791567e+00_DP
        xtab(9) =    0.632876153031869677662404854444e+00_DP
        xtab(10) =   0.819279321644006678348641581717e+00_DP
        xtab(11) =   0.944899272222882223407580138303e+00_DP
        xtab(12) =   1.0e+00_DP
      case (13)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.953309846642163911896905464755e+00_DP
        xtab(3) =  - 0.846347564651872316865925607099e+00_DP
        xtab(4) =  - 0.686188469081757426072759039566e+00_DP
        xtab(5) =  - 0.482909821091336201746937233637e+00_DP
        xtab(6) =  - 0.249286930106239992568673700374e+00_DP
        xtab(7) =    0.0e+00_DP
        xtab(8) =    0.249286930106239992568673700374e+00_DP
        xtab(9) =    0.482909821091336201746937233637e+00_DP
        xtab(10) =   0.686188469081757426072759039566e+00_DP
        xtab(11) =   0.846347564651872316865925607099e+00_DP
        xtab(12) =   0.953309846642163911896905464755e+00_DP
        xtab(13) =   1.0e+00_DP
      case (14)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.959935045267260901355100162015e+00_DP
        xtab(3) =  - 0.867801053830347251000220202908e+00_DP
        xtab(4) =  - 0.728868599091326140584672400521e+00_DP
        xtab(5) =  - 0.550639402928647055316622705859e+00_DP
        xtab(6) =  - 0.342724013342712845043903403642e+00_DP
        xtab(7) =  - 0.116331868883703867658776709736e+00_DP
        xtab(8) =    0.116331868883703867658776709736e+00_DP
        xtab(9) =    0.342724013342712845043903403642e+00_DP
        xtab(10) =   0.550639402928647055316622705859e+00_DP
        xtab(11) =   0.728868599091326140584672400521e+00_DP
        xtab(12) =   0.867801053830347251000220202908e+00_DP
        xtab(13) =   0.959935045267260901355100162015e+00_DP
        xtab(14) =   1.0e+00_DP
      case (15)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.965245926503838572795851392070e+00_DP
        xtab(3) =  - 0.885082044222976298825401631482e+00_DP
        xtab(4) =  - 0.763519689951815200704118475976e+00_DP
        xtab(5) =  - 0.606253205469845711123529938637e+00_DP
        xtab(6) =  - 0.420638054713672480921896938739e+00_DP
        xtab(7) =  - 0.215353955363794238225679446273e+00_DP
        xtab(8) =    0.0e+00_DP
        xtab(9) =    0.215353955363794238225679446273e+00_DP
        xtab(10) =   0.420638054713672480921896938739e+00_DP
        xtab(11) =   0.606253205469845711123529938637e+00_DP
        xtab(12) =   0.763519689951815200704118475976e+00_DP
        xtab(13) =   0.885082044222976298825401631482e+00_DP
        xtab(14) =   0.965245926503838572795851392070e+00_DP
        xtab(15) =   1.0e+00_DP
      case (16)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.969568046270217932952242738367e+00_DP
        xtab(3) =  - 0.899200533093472092994628261520e+00_DP
        xtab(4) =  - 0.792008291861815063931088270963e+00_DP
        xtab(5) =  - 0.652388702882493089467883219641e+00_DP
        xtab(6) =  - 0.486059421887137611781890785847e+00_DP
        xtab(7) =  - 0.299830468900763208098353454722e+00_DP
        xtab(8) =  - 0.101326273521949447843033005046e+00_DP
        xtab(9) =    0.101326273521949447843033005046e+00_DP
        xtab(10) =   0.299830468900763208098353454722e+00_DP
        xtab(11) =   0.486059421887137611781890785847e+00_DP
        xtab(12) =   0.652388702882493089467883219641e+00_DP
        xtab(13) =   0.792008291861815063931088270963e+00_DP
        xtab(14) =   0.899200533093472092994628261520e+00_DP
        xtab(15) =   0.969568046270217932952242738367e+00_DP
        xtab(16) =   1.0e+00_DP
      case (17)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.973132176631418314156979501874e+00_DP
        xtab(3) =  - 0.910879995915573595623802506398e+00_DP
        xtab(4) =  - 0.815696251221770307106750553238e+00_DP
        xtab(5) =  - 0.691028980627684705394919357372e+00_DP
        xtab(6) =  - 0.541385399330101539123733407504e+00_DP
        xtab(7) =  - 0.372174433565477041907234680735e+00_DP
        xtab(8) =  - 0.189511973518317388304263014753e+00_DP
        xtab(9) =    0.0e+00_DP
        xtab(10) =   0.189511973518317388304263014753e+00_DP
        xtab(11) =   0.372174433565477041907234680735e+00_DP
        xtab(12) =   0.541385399330101539123733407504e+00_DP
        xtab(13) =   0.691028980627684705394919357372e+00_DP
        xtab(14) =   0.815696251221770307106750553238e+00_DP
        xtab(15) =   0.910879995915573595623802506398e+00_DP
        xtab(16) =   0.973132176631418314156979501874e+00_DP
        xtab(17) =   1.0e+00_DP
      case (18)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.976105557412198542864518924342e+00_DP
        xtab(3) =  - 0.920649185347533873837854625431e+00_DP
        xtab(4) =  - 0.835593535218090213713646362328e+00_DP
        xtab(5) =  - 0.723679329283242681306210365302e+00_DP
        xtab(6) =  - 0.588504834318661761173535893194e+00_DP
        xtab(7) =  - 0.434415036912123975342287136741e+00_DP
        xtab(8) =  - 0.266362652878280984167665332026e+00_DP
        xtab(9) =  - 0.897490934846521110226450100886e-01_DP
        xtab(10) =   0.897490934846521110226450100886e-01_DP
        xtab(11) =   0.266362652878280984167665332026e+00_DP
        xtab(12) =   0.434415036912123975342287136741e+00_DP
        xtab(13) =   0.588504834318661761173535893194e+00_DP
        xtab(14) =   0.723679329283242681306210365302e+00_DP
        xtab(15) =   0.835593535218090213713646362328e+00_DP
        xtab(16) =   0.920649185347533873837854625431e+00_DP
        xtab(17) =   0.976105557412198542864518924342e+00_DP
        xtab(18) =   1.0e+00_DP
      case (19)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.978611766222080095152634063110e+00_DP
        xtab(3) =  - 0.928901528152586243717940258797e+00_DP
        xtab(4) =  - 0.852460577796646093085955970041e+00_DP
        xtab(5) =  - 0.751494202552613014163637489634e+00_DP
        xtab(6) =  - 0.628908137265220497766832306229e+00_DP
        xtab(7) =  - 0.488229285680713502777909637625e+00_DP
        xtab(8) =  - 0.333504847824498610298500103845e+00_DP
        xtab(9) =  - 0.169186023409281571375154153445e+00_DP
        xtab(10) =   0.0e+00_DP
        xtab(11) =   0.169186023409281571375154153445e+00_DP
        xtab(12) =   0.333504847824498610298500103845e+00_DP
        xtab(13) =   0.488229285680713502777909637625e+00_DP
        xtab(14) =   0.628908137265220497766832306229e+00_DP
        xtab(15) =   0.751494202552613014163637489634e+00_DP
        xtab(16) =   0.852460577796646093085955970041e+00_DP
        xtab(17) =   0.928901528152586243717940258797e+00_DP
        xtab(18) =   0.978611766222080095152634063110e+00_DP
        xtab(19) =   1.0e+00_DP
      case (20)
        xtab(1) =  - 1.0e+00_DP
        xtab(2) =  - 0.980743704893914171925446438584e+00_DP
        xtab(3) =  - 0.935934498812665435716181584931e+00_DP
        xtab(4) =  - 0.866877978089950141309847214616e+00_DP
        xtab(5) =  - 0.775368260952055870414317527595e+00_DP
        xtab(6) =  - 0.663776402290311289846403322971e+00_DP
        xtab(7) =  - 0.534992864031886261648135961829e+00_DP
        xtab(8) =  - 0.392353183713909299386474703816e+00_DP
        xtab(9) =  - 0.239551705922986495182401356927e+00_DP
        xtab(10) = - 0.805459372388218379759445181596e-01_DP
        xtab(11) =   0.805459372388218379759445181596e-01_DP
        xtab(12) =   0.239551705922986495182401356927e+00_DP
        xtab(13) =   0.392353183713909299386474703816e+00_DP
        xtab(14) =   0.534992864031886261648135961829e+00_DP
        xtab(15) =   0.663776402290311289846403322971e+00_DP
        xtab(16) =   0.775368260952055870414317527595e+00_DP
        xtab(17) =   0.866877978089950141309847214616e+00_DP
        xtab(18) =   0.935934498812665435716181584931e+00_DP
        xtab(19) =   0.980743704893914171925446438584e+00_DP
        xtab(20) =   1.0e+00_DP
      case (21)
        xtab(1) =   -1._DP
        xtab(2) =   -0.98257229660454802823448127655541_DP
        xtab(3) =   -0.94197629695974553429610265066144_DP
        xtab(4) =   -0.87929475532359046445115359630494_DP
        xtab(5) =   -0.79600192607771240474431258966036_DP
        xtab(6) =   -0.69405102606222323262731639319467_DP
        xtab(7) =   -0.57583196026183068692702187033809_DP
        xtab(8) =   -0.44411578327900210119451634960735_DP
        xtab(9) =   -0.30198985650876488727535186785875_DP
        xtab(10) =  -0.15278551580218546600635832848567_DP
        xtab(11) =   0._DP
        xtab(12) =   0.15278551580218546600635832848567_DP
        xtab(13) =   0.30198985650876488727535186785875_DP
        xtab(14) =   0.44411578327900210119451634960735_DP
        xtab(15) =   0.57583196026183068692702187033809_DP
        xtab(16) =   0.69405102606222323262731639319467_DP
        xtab(17) =   0.79600192607771240474431258966036_DP
        xtab(18) =   0.87929475532359046445115359630494_DP
        xtab(19) =   0.94197629695974553429610265066144_DP
        xtab(20) =   0.98257229660454802823448127655541_DP
        xtab(21) =   1._DP
      case (22)
        xtab(1) =   -1._DP
        xtab(2) =   -0.98415243845764617655228962221207_DP
        xtab(3) =   -0.94720428399922868052421376661573_DP
        xtab(4) =   -0.89006229019090447052965782577909_DP
        xtab(5) =   -0.81394892761192113604544184805614_DP
        xtab(6) =   -0.72048723996120215811988189639847_DP
        xtab(7) =   -0.61166943828425897122621160586993_DP
        xtab(8) =   -0.48981487518990234980875123568327_DP
        xtab(9) =   -0.35752071013891953806095728024018_DP
        xtab(10) =  -0.21760658515928504178795509346539_DP
        xtab(11) =  -0.73054540010898334761088790464107e-1_DP
        xtab(12) =   0.73054540010898334761088790464107e-1_DP
        xtab(13) =   0.21760658515928504178795509346539_DP
        xtab(14) =   0.35752071013891953806095728024018_DP
        xtab(15) =   0.48981487518990234980875123568327_DP
        xtab(16) =   0.61166943828425897122621160586993_DP
        xtab(17) =   0.72048723996120215811988189639847_DP
        xtab(18) =   0.81394892761192113604544184805614_DP
        xtab(19) =   0.89006229019090447052965782577909_DP
        xtab(20) =   0.94720428399922868052421376661573_DP
        xtab(21) =   0.98415243845764617655228962221207_DP
        xtab(22) =   1._DP
      case default
        errcode=1001
      end select
      return
      end subroutine get1Dintpolpoints
