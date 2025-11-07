module m_stl
  implicit none

  private
  public :: t_stl
  public :: STL_Init
  public :: STL_Finalize
  ! public :: STL_Add_facet
  public :: STL_Remove_duplicate_vertex
  public :: STL_Read_ascii
  public :: STL_Write_ascii
  public :: STL_Generate_cube

  type t_stl
    integer :: nfacet, nvertex
    real(kind=4), pointer :: normal(:,:) => null()
    real(kind=4), pointer :: vertex(:,:) => null()
    integer, pointer :: v_on_f(:,:) => null()
    character(len=80) :: header
  end type t_stl

  integer, parameter :: INIT_NFACET = 8

contains

  subroutine STL_Init(stl)
    type(t_stl), intent(out) :: stl
    stl%nfacet = 0
    stl%nvertex = 0
    allocate( stl%normal(3,INIT_NFACET) )
    allocate( stl%vertex(3,3*INIT_NFACET) )
    allocate( stl%v_on_f(3,INIT_NFACET) )
  end subroutine STL_Init

  subroutine STL_Finalize(stl)
    type(t_stl), intent(inout) :: stl
    if( associated( stl%normal ) ) deallocate( stl%normal )
    if( associated( stl%vertex ) ) deallocate( stl%vertex )
    if( associated( stl%v_on_f ) ) deallocate( stl%v_on_f )
    stl%nfacet = 0
    stl%nvertex = 0
  end subroutine STL_Finalize

  subroutine check_normal(normal, vertex)
    real(kind=4), intent(in) :: normal(3)
    real(kind=4), intent(in) :: vertex(3,3)
    real(kind=4) :: v12(3), v13(3), xp(3), len, dp
    v12(1:3) = vertex(1:3,2) - vertex(1:3,1)
    v13(1:3) = vertex(1:3,3) - vertex(1:3,1)
    xp(1) = v12(2) * v13(3) - v12(3) * v13(2)
    xp(2) = v12(3) * v13(1) - v12(1) * v13(3)
    xp(3) = v12(1) * v13(2) - v12(2) * v13(1)
    len = sqrt(xp(1)**2 + xp(2)**2 + xp(3)**2)
    xp(1:3) = xp(1:3) / len
    dp = dot_product(normal, xp)
    ! print*, 'DEBUG: check_normal: ', dp
    if( abs(dp - 1.0) > 1.0e-6 ) then
      if( dp < 1.0e-3 ) then
        print*, 'Warning: inconsistent surface normal', dp
      else
        print*, 'Warning: unnormalized normal', dp
      endif
    endif
  end subroutine check_normal

  subroutine STL_Add_facet(stl, normal, vertex)
    type(t_stl), intent(inout) :: stl
    real(kind=4), intent(in) :: normal(3)
    real(kind=4), intent(in) :: vertex(3,3)
    real(kind=4), pointer :: tmp_normal(:,:)
    real(kind=4), pointer :: tmp_vertex(:,:)
    integer, pointer :: tmp_v_on_f(:,:)
    integer :: i

    call check_normal(normal, vertex)

    if( stl%nfacet + 1 > ubound(stl%normal, dim=2) ) then
      allocate( tmp_normal(3, ubound(stl%normal, dim=2) * 2) )
      ! tmp_normal(1:3,1:stl%nfacet) = stl%normal(1:3,1:stl%nfacet)
      do i=1,stl%nfacet
        tmp_normal(1:3,i) = stl%normal(1:3,i)
      enddo
      deallocate( stl%normal )
      stl%normal => tmp_normal

      allocate( tmp_v_on_f(3, ubound(stl%v_on_f, dim=2) * 2) )
      ! tmp_v_on_f(1:3,1:stl%nfacet) = stl%v_on_f(1:3,1:stl%nfacet)
      do i=1,stl%nfacet
        tmp_v_on_f(1:3,i) = stl%v_on_f(1:3,i)
      enddo
      deallocate( stl%v_on_f )
      stl%v_on_f => tmp_v_on_f
    endif

    if( stl%nvertex + 3 > ubound(stl%vertex, dim=2) ) then
      allocate( tmp_vertex(3, ubound(stl%vertex, dim=2) * 2) )
      ! tmp_vertex(1:3,1:stl%nvertex) = stl%vertex(1:3,1:stl%nvertex)
      do i=1,stl%nvertex
        tmp_vertex(1:3,i) = stl%vertex(1:3,i)
      enddo
      deallocate( stl%vertex )
      stl%vertex => tmp_vertex
    endif

    stl%normal(1:3,stl%nfacet+1) = normal(1:3)
    stl%vertex(1:3,stl%nvertex+1:stl%nvertex+3) = vertex(1:3,1:3)
    do i=1,3
      stl%v_on_f(i,stl%nfacet+1) = stl%nvertex+i
    enddo

    stl%nfacet = stl%nfacet + 1
    stl%nvertex = stl%nvertex + 3
  end subroutine STL_Add_facet

  integer function compare_vertex( vtx1, vtx2 )
    real(kind=4), intent(in) :: vtx1(3), vtx2(3)
    real(kind=4) :: diff(3)
    real(kind=4), parameter :: eps = 1.e-13
    diff(1:3) = vtx1(1:3) - vtx2(1:3)
    if( diff(1) < -eps ) then
      compare_vertex = -1; return
    else if( diff(1) >  eps ) then
      compare_vertex =  1; return
    else if( diff(2) < -eps ) then
      compare_vertex = -1; return
    else if( diff(2) >  eps ) then
      compare_vertex =  1; return
    else if( diff(3) < -eps ) then
      compare_vertex = -1; return
    else if( diff(3) >  eps ) then
      compare_vertex =  1; return
    else
      compare_vertex = 0; return
    endif
  end function compare_vertex

  recursive subroutine sort_vertex( vertex, istart, iend )
    real(kind=4), intent(inout) :: vertex(:,:)
    integer, intent(in) :: istart, iend
    real(kind=4) :: pivot(3), temp(3)
    integer :: ileft, iright, icmp
    if( istart >= iend ) return
    pivot(1:3) = vertex(1:3,istart)
    ileft = istart + 1
    iright = iend
    do
      do while( ileft < iend )
        icmp = compare_vertex( vertex(1:3,ileft), pivot(1:3) )
        if( icmp >= 0 ) exit
        ileft = ileft + 1
      enddo
      do while( istart < iright )
        icmp = compare_vertex( vertex(1:3,iright), pivot(1:3) )
        if( icmp < 0 ) exit
        iright = iright - 1
      enddo
      if( ileft > iright ) exit
      temp(1:3) = vertex(1:3,ileft)
      vertex(1:3,ileft) = vertex(1:3,iright)
      vertex(1:3,iright) = temp(1:3)
      ileft = ileft + 1
      iright = iright - 1
      if( ileft > iright ) exit
    enddo
    ileft = ileft - 1
    temp(1:3) = vertex(1:3,ileft)
    vertex(1:3,ileft) = vertex(1:3,istart)
    vertex(1:3,istart) = temp(1:3)
    call sort_vertex( vertex, istart, ileft-1 )
    call sort_vertex( vertex, iright+1, iend )
  end subroutine sort_vertex

  subroutine uniq_vertex( vertex, nvertex, nvertex_new )
    real(kind=4), intent(inout) :: vertex(:,:)
    integer, intent(in) :: nvertex
    integer, intent(out) :: nvertex_new
    integer :: ndup, i, icmp
    ndup = 0
    do i=2,nvertex
      icmp = compare_vertex( vertex(1:3,i), vertex(1:3,i-1-ndup) )
      if( icmp == 0 ) then
        ndup = ndup + 1
      else if( ndup > 0 ) then
        vertex(1:3,i-ndup) = vertex(1:3,i)
      endif
    enddo
    nvertex_new = nvertex - ndup
  end subroutine uniq_vertex

  recursive subroutine bsearch_vertex( vertex, istart, iend, vtx, id )
    real(kind=4), intent(in) :: vertex(:,:)
    integer, intent(in) :: istart, iend
    real(kind=4), intent(in) :: vtx(3)
    integer, intent(out) :: id
    integer :: icenter, icmp
    if( istart > iend ) then
      id = -1
      return
    endif
    icenter = (istart + iend) / 2
    icmp = compare_vertex( vtx, vertex(1:3,icenter) )
    if( icmp < 0 ) then
      call bsearch_vertex( vertex, istart, icenter-1, vtx, id )
    else if( icmp > 0 ) then
      call bsearch_vertex( vertex, icenter+1, iend, vtx, id )
    else
      id = icenter
    endif
  end subroutine bsearch_vertex

  subroutine STL_Remove_duplicate_vertex(stl)
    type(t_stl), intent(inout) :: stl

    real(kind=4), pointer :: tmp_vertex(:,:)
    integer :: nvertex_new, i, j, id

    allocate( tmp_vertex(3,1:stl%nvertex) )
    ! tmp_vertex(1:3,1:stl%nvertex) = stl%vertex(1:3,1:stl%nvertex)
    do i=1,stl%nvertex
      tmp_vertex(1:3,i) = stl%vertex(1:3,i)
    enddo

    call sort_vertex( tmp_vertex, 1, stl%nvertex )
    call uniq_vertex( tmp_vertex, stl%nvertex, nvertex_new )
    ! write(0,*) 'DEBUG: nvertex,nvertex_new:',stl%nvertex,nvertex_new

    do i=1,stl%nfacet
      do j=1,3
        call bsearch_vertex( tmp_vertex, 1, nvertex_new, stl%vertex(1:3,stl%v_on_f(j,i)), id )
        if( id < 0 ) then
          print*, 'ERROR: vertex not found'; stop
        endif
        stl%v_on_f(j,i) = id
      enddo
    enddo

    deallocate( stl%vertex )
    stl%vertex => tmp_vertex
    stl%nvertex = nvertex_new
  end subroutine STL_Remove_duplicate_vertex

  subroutine STL_Read_ascii(stl, filename, ierr)
    type(t_stl), intent(out) :: stl
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ierr

    character(len=256) :: record
    integer :: iunit, nline, i
    real(kind=4) :: normal(3), vertex(3,3)

    ierr = 0

    call STL_Finalize(stl)
    call STL_Init(stl)

    open( newunit=iunit, file=filename, status='old', action='read', iostat=ierr )
    if( ierr /= 0 ) then
      print*, trim(filename), ' does not exist'; stop
    endif
    nline = 0

    read( iunit, fmt='(a)', iostat=ierr ) record
    if( ierr /= 0 ) return
    nline = nline + 1
    record = adjustl(record)
    if( record(1:5) /= 'solid' ) then
      print*, 'Input file is not in STL ascii format'
      ierr = 1
      return
    endif
    stl%header = ''
    if( record(6:6) == ' ' ) then
      stl%header(1:80) = adjustl(record(7:86))
    endif

    read_facet: do
      read( iunit, fmt='(a)', iostat=ierr ) record
      if( ierr /= 0 ) then
        print*, 'File truncated at ', trim(filename), ':', nline
        ierr = 1
        return
      endif
      nline = nline + 1
      record = adjustl(record)
      if( record(1:13) /= 'facet normal ' ) then
        exit read_facet
      endif
      read( record(14:), fmt=*, iostat=ierr ) normal(1:3)
      if( ierr /= 0 ) then
        print*, 'Format error in ', trim(filename), ':', nline
        return
      endif

      read( iunit, fmt='(a)', iostat=ierr ) record
      if( ierr /= 0 ) then
        print*, 'File truncated at ', trim(filename), ':', nline
        ierr = 1
        return
      endif
      nline = nline + 1
      record = adjustl(record)
      if( record(1:10) /= 'outer loop' ) then
        print*, 'Format error in ', trim(filename), ':', nline
        ierr = 1
        return
      endif

      do i=1,3
        read( iunit, fmt='(a)', iostat=ierr ) record
        if( ierr /= 0 ) then
          print*, 'File truncated at ', trim(filename), ':', nline
          ierr = 1
          return
        endif
        nline = nline + 1
        record = adjustl(record)
        if( record(1:7) /= 'vertex ' ) then
          print*, 'Format error in ', trim(filename), ':', nline
          ierr = 1
          return
        endif
        read( record(8:), fmt=*, iostat=ierr ) vertex(1:3,i)
        if( ierr /= 0 ) then
          print*, 'Format error in ', trim(filename), ':', nline
          return
        endif
      enddo

      read( iunit, fmt='(a)', iostat=ierr ) record
      if( ierr /= 0 ) then
        print*, 'File truncated at ', trim(filename), ':', nline
        ierr = 1
        return
      endif
      nline = nline + 1
      record = adjustl(record)
      if( record(1:7) /= 'endloop' ) then
        print*, 'Format error in ', trim(filename), ':', nline
        ierr = 1
        return
      endif

      read( iunit, fmt='(a)', iostat=ierr ) record
      if( ierr /= 0 ) then
        print*, 'File truncated at ', trim(filename), ':', nline
        ierr = 1
        return
      endif
      nline = nline + 1
      record = adjustl(record)
      if( record(1:8) /= 'endfacet' ) then
        print*, 'Format error in ', trim(filename), ':', nline
        ierr = 1
        return
      endif

      call STL_Add_facet(stl, normal, vertex)
    enddo read_facet

    if( record(1:8) /= 'endsolid' ) then
      print*, 'Format error in ', trim(filename), ':', nline
      ierr = 1
      return
    endif

    close( iunit )

  end subroutine STL_Read_ascii

  subroutine STL_Write_ascii(stl, filename, ierr)
    type(t_stl), intent(in) :: stl
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ierr

    integer :: iunit, i, j

    open( newunit=iunit, file=filename, status='replace', action='write', iostat=ierr )
    if( ierr /= 0 ) then
      print*, 'Failed to opne file ', trim(filename), ' for write'; stop
    endif

    write( iunit, '(a,a)' ) 'solid ',trim(stl%header)

    do i=1,stl%nfacet
      write( iunit, '(a,1p3e14.6)' ) '  facet normal ', stl%normal(1:3,i)
      write( iunit, '(a)' ) '    outer loop'
      do j=1,3
        write( iunit, '(a,1p3e14.6)' ) '      vertex ', stl%vertex(1:3,stl%v_on_f(j,i))
      enddo
      write( iunit, '(a)' ) '    endloop'
      write( iunit, '(a)' ) '  endfacet'
    enddo

    write( iunit, '(a)' ) 'endsolid'

    close( iunit )

  end subroutine STL_Write_ascii

  subroutine STL_Generate_cube(stl, edgelen, ndiv)
    type(t_stl), intent(out) :: stl
    real(kind=4), intent(in) :: edgelen
    integer, intent(in) :: ndiv

    real(kind=4) :: normal(3), vertex(3,3)
    integer :: i, j
    real(kind=4) :: x0, x1, x2, y0, y1, y2, z0, z1, z2

    call STL_Finalize(stl)
    call STL_Init(stl)

    write(stl%header, '(a,1p1e13.6,a,i0)') 'cube with edgelen=',edgelen,', ndiv=',ndiv

    x0 = 0; y0 = 0; z0 = 0
    ! x0 = -edgelen/2; y0 = -edgelen/2; z0 = -edgelen/2

    ! x = x0
    vertex(1,1:3) = x0
    normal(1) = -1; normal(2) = 0; normal(3) = 0
    do j=1,ndiv
      z1 = z0 + edgelen * (j-1) / ndiv
      z2 = z0 + edgelen * j / ndiv
      do i=1,ndiv
        y1 = y0 + edgelen * (i-1) / ndiv
        y2 = y0 + edgelen * i / ndiv
        vertex(2,1) = y1; vertex(3,1) = z1
        vertex(2,2) = y1; vertex(3,2) = z2
        vertex(2,3) = y2; vertex(3,3) = z1
        call STL_Add_facet(stl, normal, vertex)
        vertex(1:3,1) = vertex(1:3,3)
        vertex(2,3) = y2; vertex(3,3) = z2
        call STL_Add_facet(stl, normal, vertex)
      enddo
    enddo

    ! x = x0 + edgelen
    vertex(1,1:3) = x0 + edgelen
    normal(1) = 1; normal(2) = 0; normal(3) = 0
    do j=1,ndiv
      z1 = z0 + edgelen * (j-1) / ndiv
      z2 = z0 + edgelen * j / ndiv
      do i=1,ndiv
        y1 = y0 + edgelen * (i-1) / ndiv
        y2 = y0 + edgelen * i / ndiv
        vertex(2,1) = y1; vertex(3,1) = z1
        vertex(2,2) = y2; vertex(3,2) = z1
        vertex(2,3) = y1; vertex(3,3) = z2
        call STL_Add_facet(stl, normal, vertex)
        vertex(1:3,1) = vertex(1:3,2)
        vertex(2,2) = y2; vertex(3,2) = z2
        call STL_Add_facet(stl, normal, vertex)
      enddo
    enddo

    ! y = y0
    vertex(2,1:3) = y0
    normal(1) = 0; normal(2) = -1; normal(3) = 0
    do j=1,ndiv
      x1 = x0 + edgelen * (j-1) / ndiv
      x2 = x0 + edgelen * j / ndiv
      do i=1,ndiv
        z1 = z0 + edgelen * (i-1) / ndiv
        z2 = z0 + edgelen * i / ndiv
        vertex(1,1) = x1; vertex(3,1) = z1
        vertex(1,2) = x2; vertex(3,2) = z1
        vertex(1,3) = x1; vertex(3,3) = z2
        call STL_Add_facet(stl, normal, vertex)
        vertex(1:3,1) = vertex(1:3,3)
        vertex(1,3) = x2; vertex(3,3) = z2
        call STL_Add_facet(stl, normal, vertex)
      enddo
    enddo

    ! y = y0 + edgelen
    vertex(2,1:3) = y0 + edgelen
    normal(1) = 0; normal(2) = 1; normal(3) = 0
    do j=1,ndiv
      x1 = x0 + edgelen * (j-1) / ndiv
      x2 = x0 + edgelen * j / ndiv
      do i=1,ndiv
        z1 = z0 + edgelen * (i-1) / ndiv
        z2 = z0 + edgelen * i / ndiv
        vertex(1,1) = x1; vertex(3,1) = z1
        vertex(1,2) = x1; vertex(3,2) = z2
        vertex(1,3) = x2; vertex(3,3) = z1
        call STL_Add_facet(stl, normal, vertex)
        vertex(1:3,1) = vertex(1:3,2)
        vertex(1,2) = x2; vertex(3,2) = z2
        call STL_Add_facet(stl, normal, vertex)
      enddo
    enddo

    ! z = z0
    vertex(3,1:3) = z0
    normal(1) = 0; normal(2) = 0; normal(3) = -1
    do j=1,ndiv
      y1 = y0 + edgelen * (j-1) / ndiv
      y2 = y0 + edgelen * j / ndiv
      do i=1,ndiv
        x1 = x0 + edgelen * (i-1) / ndiv
        x2 = x0 + edgelen * i / ndiv
        vertex(1,1) = x1; vertex(2,1) = y1
        vertex(1,2) = x1; vertex(2,2) = y2
        vertex(1,3) = x2; vertex(2,3) = y1
        call STL_Add_facet(stl, normal, vertex)
        vertex(1:3,1) = vertex(1:3,3)
        vertex(1,3) = x2; vertex(2,3) = y2
        call STL_Add_facet(stl, normal, vertex)
      enddo
    enddo

    ! z = z0 + edgelen
    vertex(3,1:3) = z0 + edgelen
    normal(1) = 0; normal(2) = 0; normal(3) = 1
    do j=1,ndiv
      y1 = y0 + edgelen * (j-1) / ndiv
      y2 = y0 + edgelen * j / ndiv
      do i=1,ndiv
        x1 = x0 + edgelen * (i-1) / ndiv
        x2 = x0 + edgelen * i / ndiv
        vertex(1,1) = x1; vertex(2,1) = y1
        vertex(1,2) = x2; vertex(2,2) = y1
        vertex(1,3) = x1; vertex(2,3) = y2
        call STL_Add_facet(stl, normal, vertex)
        vertex(1:3,1) = vertex(1:3,2)
        vertex(1,2) = x2; vertex(2,2) = y2
        call STL_Add_facet(stl, normal, vertex)
      enddo
    enddo

  end subroutine STL_Generate_cube

end module m_stl
