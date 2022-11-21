      subroutine rbtxt(txtgeb,anzgeb,alrb,btrb,donetx,lplayt,layrb,nnat)
      use feminterface, only: strtok
      use femtypes
      implicit none
      integer (I4B) :: anzgeb, lplayt(:), nnat
      integer (I4B) :: layrb(:,:), donetx(:)
      complex (DPC) :: alrb(:,:), btrb(:,:)
      character (len=*) :: txtgeb(:)
      intent (in) :: txtgeb, anzgeb, lplayt, nnat
      intent (out) ::  alrb, btrb, donetx, layrb
!
!    $Revision: 1.7 $
!    $Date: 2011/08/10 13:16:52 $
!    $Author: juryanatzki $
!
!  Test whether a text string constitute a boundary condition.
!  For each nature, the text to decribe a boundary condition consits of up to 6 token,
!  which are sparated by a delimiter.
!  valid descriptors of boundary conditions are:
!    alpha= 1.             ;    alph = (.001,2.); beta=5.;    a= 1,1 b=2,2
!  imaginary part may be omitted, if present, the real part must be present as well
!  if beta (or b) is not present the BC is treated as a Dirichlet BC
!
!  Input:
!    txtgeb    the string of TEXT and MTEXT entities
!    anzgeb    number of texts read from the DXF-file (number of regions)
!    lplayt    layer of the text strings read from the DXF-file
!    nnat      number of natures
!
!  Output:
!    alrb      first value of boundary value descriptor
!    btrb      second value of boundary value descriptor for Neuman and general BCs
!    donetx    a flag which is 1 if the text has been recognized as a boundary condition, otherwise 0
!    layrb     type of boundary condition for branches on that layer (and for each nature)
!              0    = Dirichlet
!              1-99 = User Dirichlet
!              200  = general Neumann
!              300  = innner contour
!              400  = non-visble contour
!              1000 = comment layer
!              1001 = not recognized layer, treated as a comment layer
!
!  local variables
      integer (I4B) :: i, j, ios, idx, lay
      real (DP) :: ar, ai, br, bi
      character (len=1000) :: string
      character (len=80) :: token
      character (len=20) :: delim
!
!  Delimiter for tokenizing a string
!      delim = '=,;:()?'//char(9)//char(34)//char(13)  ! char(9)= horizontal tab; char(13)= CR; char(34)="
      delim=' =,;:"{}()!@#$%^&*'
!
      do i=1,anzgeb
        donetx(i) = 0
!  ignore text on layers which are just comment or contour
!#######  'CONTO','GEOME' 'INVIS','KOMME','COMME','BEMAS'
        string = txtgeb(i)
!  get position of first '=' sign
        idx=index(string,'=')
!  no '=' sign ?
        if (idx .eq. 0) cycle
!  get the first token 1. -> a text starting with 'a'
        lay=lplayt(i)
        call strtok(string, delim, token)
        do j=1, nnat
          if (token(1:1) .ne. 'a' .and. token(1:1) .ne. 'A' ) then
            if (j.ne.1) then
              print*,'*ERRROR: incomplete descriptor of boundary condition'
              write(*,'(a,i3,a)')' *        expecting descriptor for',nnat,' natures'
              print*,'  : ',trim(txtgeb(i))
            end if
            exit
          end if
!  get the next token 2.  -> a number, real part of A
          call strtok(char(0), delim, token)
          if (token(1:1) .ne. char(0)) then
!  convert the string:"token" to a number
            read(token,*,IOSTAT=ios) ar
            if (ios .ne. 0) then
              if (nnat .ne. 1) then
                print*,'*ERRROR: incomplete descriptor of boundary condition'
                print*,'  : ',trim(txtgeb(i))
              end if
              exit
            end if
          else
            print*,'*ERRROR: incomplete descriptor of boundary condition'
            print*,'  : ',trim(txtgeb(i))
            exit
          end if
          donetx(i) = 1
!  get the next token 3.  -> a number or a text
          call strtok(char(0), delim, token)
          if (token(1:1) .ne. char(0)) then
!  convert the string:"token" to a number
            read(token,*,IOSTAT=ios) ai
            if (ios .eq. 0) then
              alrb(lay,j) = cmplx(ar,ai)
!  get the next token 4.  -> a text starting with 'B'
              call strtok(char(0), delim, token)
            else
!  it is not a number, imaginary part of A is missing
              alrb(lay,j) = cmplx(ar,0._DP)
            end if
          else
            if (j.ne.nnat) then
              print*,'*ERRROR: incomplete descriptor of boundary condition'
              write(*,'(a,i3,a)')' *        expecting descriptor for',nnat,' natures'
              print*,'  : ',trim(txtgeb(i))
            else
!  end of string reached "B=.." is not present
              alrb(lay,j) = cmplx(ar,0._DP)
              btrb(lay,j) = cmplx(0._DP,0._DP)
              layrb(lay,j)=0
            end if
            exit
          end if
!
          if (token(1:1) .eq. 'b' .or. token(1:1) .eq. 'B' ) then
!  "B..= " is  present, it is a Neumann or general boundary condition
            layrb(lay,j)=200
!  get the next token 5.  -> a number, real part of B
            call strtok(char(0), delim, token)
            if (token(1:1) .ne. char(0)) then
!  convert the string:"token" to a number
              read(token,*,IOSTAT=ios) br
              if (ios .ne. 0) then
                print*,'*ERRROR: incomplete descriptor of boundary condition'
                print*,'  : ',trim(txtgeb(i))
                btrb(lay,j) = cmplx(0._DP,0._DP)
                exit
              end if
            else
              print*,'*ERRROR: incomplete descriptor of boundary condition'
              print*,'  : ',trim(txtgeb(i))
              exit
            end if
!  get the next token 6.  -> a number, imaginary part of B
            call strtok(char(0), delim, token)
!  convert the string:"token" to a number
            if (token(1:1) .ne. char(0)) then
              read(token,*,IOSTAT=ios) bi
              if (ios .eq. 0) then
                btrb(lay,j) = cmplx(br,bi)
!  get the first token of next nature
                if (j.ne.nnat) call strtok(char(0), delim, token)
              else
!  it is not a number, imaginary part of B is missing
                btrb(lay,j) = cmplx(br,0._DP)
              end if
            else
              if (j.ne.nnat) then
                print*,'*ERRROR: incomplete descriptor of boundary condition'
                write(*,'(a,i3,a)')' *        expecting descriptor for',nnat,' natures'
                print*,'  : ',trim(txtgeb(i))
              else
!  end of string reached imaginary part of b is not present
                btrb(lay,j) = cmplx(br,0._DP)
              end if
              exit
            end if
          else
!  "B..= " is not present, it is a Dirichlet boundary condition
            layrb(lay,j)=0
            btrb(lay,j) = cmplx(0._DP,0._DP)
          end if
        end do
      end do
!  testing
      do i=2,anzgeb
        if (donetx(i) .eq. 1) then
          do j=1,i-1
            if (donetx(j) .eq. 1) then
              if (lplayt(i) .eq. lplayt(j) ) then
                print*,'*WARNING: on the same layer multiple setting ', &
     &            'for the boundary condition was found:'
                print*,'1: ',trim(txtgeb(i))
                print*,'2: ',trim(txtgeb(j))
              end if
            end if
          end do
        end if
      end do
      return
      end subroutine rbtxt
