module m_HACApK_base_c
  use m_HACApK_base
  use iso_c_binding
  implicit none

  private
  public :: HACApK_generate_frame_blrleaf_c

  integer(c_int),pointer :: lod(:),lsp(:),lnp(:),lthr(:),lpmd(:)
  real(c_double),pointer :: param(:),time(:)

  interface
    function cHACApK_copy_lcontrol_init() bind(c,name='cHACApK_copy_lcontrol_init')
      use iso_c_binding
      integer(c_int) :: cHACApK_copy_lcontrol_init
    end function cHACApK_copy_lcontrol_init

    function cHACApK_copy_lcontrol_finalize() bind(c,name='cHACApK_copy_lcontrol_finalize')
      use iso_c_binding
      integer(c_int) :: cHACApK_copy_lcontrol_finalize
    end function cHACApK_copy_lcontrol_finalize

    function cHACApK_copy_lcontrol_f2c(lod,lsp,lnp,lthr,lpmd,param,time,lf_umpi) &
        bind(c,name='cHACApK_copy_lcontrol_f2c')
      use iso_c_binding
      integer(c_int) :: cHACApK_copy_lcontrol_f2c
      type(c_ptr),value :: lod,lsp,lnp,lthr,lpmd
      type(c_ptr),value :: param,time
      integer(c_int),value :: lf_umpi
    end function cHACApK_copy_lcontrol_f2c

    function cHACApK_copy_lcontrol_c2f(lod,lsp,lnp,lthr,lpmd,param,time,lf_umpi) &
        bind(c,name='cHACApK_copy_lcontrol_c2f')
      use iso_c_binding
      integer(c_int) :: cHACApK_copy_lcontrol_c2f
      type(c_ptr) :: lod,lsp,lnp,lthr,lpmd
      type(c_ptr) :: param,time
      integer(c_int) :: lf_umpi
    end function cHACApK_copy_lcontrol_c2f

    function cHACApK_copy_leafmtxp_init() bind(c,name='cHACApK_copy_leafmtxp_init')
      use iso_c_binding
      integer(c_int) :: cHACApK_copy_leafmtxp_init
    end function cHACApK_copy_leafmtxp_init

    function cHACApK_copy_leafmtxp_finalize() bind(c,name='cHACApK_copy_leafmtxp_finalize')
      use iso_c_binding
      integer(c_int) :: cHACApK_copy_leafmtxp_finalize
    end function cHACApK_copy_leafmtxp_finalize

    function cHACApK_copy_leafmtxp_c2f_scalars(nd,nlf,nlfkt,ktmax,nbl,nlfalt,nlfl,nlft,ndlfs,ndtfs,&
        nrank_l,nrank_t,npgl) bind(c,name='cHACApK_copy_leafmtxp_c2f_scalars')
      use iso_c_binding
      integer(c_int) :: cHACApK_copy_leafmtxp_c2f_scalars
      integer(c_int) :: nd
      integer(c_int) :: nlf
      integer(c_int) :: nlfkt
      integer(c_int) :: ktmax
      integer(c_int) :: nbl
      integer(c_int) :: nlfalt
      integer(c_int) :: nlfl,nlft
      integer(c_int) :: ndlfs,ndtfs
      integer(c_int) :: nrank_l,nrank_t,npgl
    end function cHACApK_copy_leafmtxp_c2f_scalars

    function cHACApK_copy_leafmtxp_c2f_arrays(lnlfl2g,lbstrtl,lbstrtt,lbndl,lbndt,lbndlfs,lbndtfs,lbl2t) &
        bind(c,name='cHACApK_copy_leafmtxp_c2f_arrays')
      use iso_c_binding
      integer(c_int) :: cHACApK_copy_leafmtxp_c2f_arrays
      integer(c_int) :: lnlfl2g(*)
      type(c_ptr) :: lbstrtl,lbstrtt
      type(c_ptr) :: lbndl,lbndt
      type(c_ptr) :: lbndlfs,lbndtfs
      type(c_ptr) :: lbl2t
    end function cHACApK_copy_leafmtxp_c2f_arrays

    function cHACApK_copy_leafmtx_c2f_scalars(ndpth,idx,ltmtx,kt,nstrtl,ndl,nstrtt,ndt,nlf) &
        bind(c,name='cHACApK_copy_leafmtx_c2f_scalars')
      use iso_c_binding
      integer(c_int) :: cHACApK_copy_leafmtx_c2f_scalars
      integer(c_int),value :: ndpth
      integer(c_int) :: idx(*)
      integer(c_int) :: ltmtx
      integer(c_int) :: kt
      integer(c_int) :: nstrtl,ndl
      integer(c_int) :: nstrtt,ndt
      integer(c_int) :: nlf
    end function cHACApK_copy_leafmtx_c2f_scalars

    function cHACApK_generate_frame_blrleaf_if(i_bemv,gmid,lnmtx,nofc,nffc,ndim) &
        bind(c,name='cHACApK_generate_frame_blrleaf_if')
      use iso_c_binding
      integer(c_int) :: cHACApK_generate_frame_blrleaf_c
      integer(c_int),value :: i_bemv
      real(c_double) :: gmid(*)
      integer(c_int) :: lnmtx(*)
      integer(c_int),value :: nofc,nffc,ndim
    end function cHACApK_generate_frame_blrleaf_if

  end interface
contains

!***HACApK_generate_frame_blrleaf_c
  subroutine HACApK_generate_frame_blrleaf_c(st_leafmtxp,st_bemv,st_ctl,gmid,lnmtx,nofc,nffc,ndim)
    use m_HACApK_calc_entry_record
    use iso_c_binding
    !$use omp_lib
    type(st_HACApK_leafmtxp), intent(out) :: st_leafmtxp
    type(st_HACApK_calc_entry), intent(in) :: st_bemv
    type(st_HACApK_lcontrol), intent(inout) :: st_ctl
    real*8, intent(in) :: gmid(nofc,ndim)
    integer*4, intent(inout) :: lnmtx(4)
    integer*4, intent(in) :: nofc,nffc,ndim
    integer*4 :: i_bemv,ierr,il
    real(c_double) :: gmid_c(nofc*ndim)
    integer(c_int) :: lnmtx_c(0:4)

    call HACApK_calc_entry_record_set(st_bemv,i_bemv)
    ierr=cHACApK_copy_lcontrol_init()
    ierr=cHACApK_copy_leafmtxp_init()

    call HACApK_copy_lcontrol_f2c(st_ctl)
    do il=1,ndim
      gmid_c(nofc*(il-1)+1:nofc*il)=gmid(1:nofc,il)
    enddo
    lnmtx_c(1:4)=lnmtx(1:4)

    ierr=cHACApK_generate_frame_blrleaf_if(i_bemv,gmid_c,lnmtx_c,nofc,nffc,ndim)

    call HACApK_copy_leafmtxp_c2f(st_leafmtxp)
    call HACApK_copy_lcontrol_c2f(st_ctl)
    lnmtx(1:4)=lnmtx_c(1:4)

    ierr=cHACApK_copy_leafmtxp_finalize()
    ierr=cHACApK_copy_lcontrol_finalize()
    call HACApK_calc_entry_record_clear(i_bemv)
  end subroutine HACApK_generate_frame_blrleaf_c

!***HACApK_copy_lcontrol_f2c
  subroutine HACApK_copy_lcontrol_f2c(st_ctl)
    use iso_c_binding
    type(st_HACApK_lcontrol) :: st_ctl
    integer(c_int) :: nd,nrank,nthr
    ! integer(c_int),pointer :: lod(:),lsp(:),lnp(:),lthr(:),lpmd(:)
    ! real(c_double),pointer :: param(:),time(:)
    integer(c_int) :: lf_umpi
    integer*4 :: ierr

    nd=ubound(st_ctl%lod,dim=1)
    nrank=ubound(st_ctl%lsp,dim=1)
    nthr=ubound(st_ctl%lthr,dim=1)-1

    ! arrays are allocated in Fortran and used in C
    allocate(lod(0:nd),lsp(0:nrank),lnp(0:nrank),lthr(0:nthr+1),lpmd(0:50))
    allocate(param(0:100),time(0:10))

    lod(1:nd)=st_ctl%lod(1:nd)
    lsp(1:nrank)=st_ctl%lsp(1:nrank)
    lnp(1:nrank)=st_ctl%lnp(1:nrank)
    lthr(1:nthr+1)=st_ctl%lthr(1:nthr+1)
    lpmd(1:50)=st_ctl%lpmd(1:50)
    param(1:100)=st_ctl%param(1:100)
    time(1:10)=st_ctl%time(1:10)
    lf_umpi=st_ctl%lf_umpi

    ierr=cHACApK_copy_lcontrol_f2c(c_loc(lod),c_loc(lsp),c_loc(lnp),c_loc(lthr),c_loc(lpmd),&
        c_loc(param),c_loc(time),lf_umpi)

  end subroutine HACApK_copy_lcontrol_f2c

!***HACApK_copy_lcontrol_c2f
  subroutine HACApK_copy_lcontrol_c2f(st_ctl)
    use iso_c_binding
    type(st_HACApK_lcontrol) :: st_ctl
    integer(c_int) :: nd,nrank,nthr
    type(c_ptr) :: lod_p,lsp_p,lnp_p,lthr_p,lpmd_p
    type(c_ptr) :: param_p,time_p
    integer(c_int) :: lf_umpi_c
    ! integer(c_int),pointer :: lod(:),lsp(:),lnp(:),lthr(:),lpmd(:)
    ! real(c_double),pointer :: param(:),time(:)
    integer*4 :: ierr

    nd=ubound(st_ctl%lod,dim=1)
    nrank=ubound(st_ctl%lsp,dim=1)
    nthr=ubound(st_ctl%lthr,dim=1)-1

    ierr=cHACApK_copy_lcontrol_c2f(lod_p,lsp_p,lnp_p,lthr_p,lpmd_p,param_p,time_p,lf_umpi_c)

    ! call c_f_pointer(lod_p,lod,[nd+1])
    ! call c_f_pointer(lsp_p,lsp,[nrank+1])
    ! call c_f_pointer(lnp_p,lnp,[nrank+1])
    ! call c_f_pointer(lthr_p,lthr,[nthr+2])
    ! call c_f_pointer(lpmd_p,lpmd,[51])
    ! call c_f_pointer(param_p,param,[101])
    ! call c_f_pointer(time_p,time,[11])
    ! st_ctl%lod(1:nd)=lod(2:nd+1)
    ! st_ctl%lsp(1:nrank)=lsp(2:nrank+1)
    ! st_ctl%lnp(1:nrank)=lnp(2:nrank+1)
    ! st_ctl%lthr(1:nthr+1)=lthr(2:nthr+2)
    ! st_ctl%lpmd(1:50)=lpmd(2:51)
    ! st_ctl%param(1:100)=param(2:101)
    ! st_ctl%time(1:10)=time(2:11)

    st_ctl%lod(1:nd)=lod(1:nd)
    st_ctl%lsp(1:nrank)=lsp(1:nrank)
    st_ctl%lnp(1:nrank)=lnp(1:nrank)
    st_ctl%lthr(1:nthr+1)=lthr(1:nthr+1)
    st_ctl%lpmd(1:50)=lpmd(1:50)
    st_ctl%param(1:100)=param(1:100)
    st_ctl%time(1:10)=time(1:10)

    st_ctl%lf_umpi=lf_umpi_c

    ! arrays used in C are allocated in Fortran; deallocate them here
    deallocate(lod,lsp,lnp,lthr,lpmd)
    deallocate(param,time)
  end subroutine HACApK_copy_lcontrol_c2f

!***HACApK_copy_leafmtxp_c2f
  subroutine HACApK_copy_leafmtxp_c2f(st_leafmtxp)
    use iso_c_binding
    type(st_HACApK_leafmtxp) :: st_leafmtxp
    integer(c_int) :: nd_c,nlf_c,nlfkt_c,ktmax_c,nbl_c,nlfalt_c,nlfl_c,nlft_c,ndlfs_c,ndtfs_c
    integer(c_int) :: nrank_l_c,nrank_t_c,npgl_c
    type(c_ptr) :: lbstrtl_p,lbstrtt_p,lbndl_p,lbndt_p,lbndlfs_p,lbndtfs_p,lbl2t_p
    integer(c_int),allocatable :: lnlfl2g_c(:)
    integer(c_int),pointer :: lbstrtl_c(:),lbstrtt_c(:)
    integer(c_int),pointer :: lbndl_c(:),lbndt_c(:)
    integer(c_int),pointer :: lbndlfs_c(:),lbndtfs_c(:)
    integer(c_int),pointer :: lbl2t_c(:)
    integer*4 :: ierr,il

    ierr=cHACApK_copy_leafmtxp_c2f_scalars(nd_c,nlf_c,nlfkt_c,ktmax_c,nbl_c,nlfalt_c,nlfl_c,nlft_c,ndlfs_c,ndtfs_c,&
        nrank_l_c,nrank_t_c,npgl_c)
    st_leafmtxp%nd=nd_c
    st_leafmtxp%nlf=nlf_c
    st_leafmtxp%nlfkt=nlfkt_c
    st_leafmtxp%ktmax=ktmax_c
    st_leafmtxp%nbl=nbl_c
    st_leafmtxp%nlfalt=nlfalt_c
    st_leafmtxp%nlfl=nlfl_c
    st_leafmtxp%nlft=nlft_c
    st_leafmtxp%ndlfs=ndlfs_c
    st_leafmtxp%ndtfs=ndtfs_c

    allocate(lnlfl2g_c(nlft_c*nlfl_c))

    ierr=cHACApK_copy_leafmtxp_c2f_arrays(lnlfl2g_c,lbstrtl_p,lbstrtt_p,lbndl_p,lbndt_p,lbndlfs_p,lbndtfs_p,lbl2t_p)

    call c_f_pointer(lbstrtl_p,lbstrtl_c,[nlfalt_c+2])
    call c_f_pointer(lbstrtt_p,lbstrtt_c,[nlfalt_c+2])
    call c_f_pointer(lbndl_p,lbndl_c,[nlfalt_c+1])
    call c_f_pointer(lbndt_p,lbndt_c,[nlfalt_c+1])
    call c_f_pointer(lbndlfs_p,lbndlfs_c,[nrank_l_c])
    call c_f_pointer(lbndtfs_p,lbndtfs_c,[nrank_t_c])
    call c_f_pointer(lbl2t_p,lbl2t_c,[npgl_c])

    allocate(st_leafmtxp%lnlfl2g(nlft_c,nlfl_c))
    allocate(st_leafmtxp%lbstrtl(nlfalt_c+1),st_leafmtxp%lbstrtt(nlfalt_c+1))
    allocate(st_leafmtxp%lbndl(nlfalt_c),st_leafmtxp%lbndt(nlfalt_c))
    allocate(st_leafmtxp%lbndlfs(0:nrank_l_c-1),st_leafmtxp%lbndtfs(0:nrank_t_c-1))
    allocate(st_leafmtxp%lbl2t(0:npgl_c-1))

    do il=1,nlfl_c
      st_leafmtxp%lnlfl2g(1:nlft_c,il)=lnlfl2g_c(nlft_c*(il-1)+1:nlft_c*il)
    enddo
    st_leafmtxp%lbstrtl(1:nlfalt_c+1)=lbstrtl_c(2:nlfalt_c+2)
    st_leafmtxp%lbstrtt(1:nlfalt_c+1)=lbstrtt_c(2:nlfalt_c+2)
    st_leafmtxp%lbndl(1:nlfalt_c)=lbndl_c(2:nlfalt_c+1)
    st_leafmtxp%lbndt(1:nlfalt_c)=lbndt_c(2:nlfalt_c+1)
    st_leafmtxp%lbndlfs(0:nrank_l_c-1)=lbndlfs_c(1:nrank_l_c)
    st_leafmtxp%lbndtfs(0:nrank_t_c-1)=lbndtfs_c(1:nrank_t_c)
    st_leafmtxp%lbl2t(0:npgl_c-1)=lbl2t_c(1:npgl_c)

    deallocate(lnlfl2g_c)

    call HACApK_copy_leafmtx_c2f(st_leafmtxp%st_lf,st_leafmtxp%nlf,0)
  end subroutine HACApK_copy_leafmtxp_c2f

!***HACApK_copy_leafmtx_c2f
  recursive subroutine HACApK_copy_leafmtx_c2f(st_lf,nlf,ndpth,idx)
    use iso_c_binding
    type(st_HACApK_leafmtx),pointer :: st_lf(:)
    integer*4, intent(in) :: nlf
    integer(c_int),intent(in) :: ndpth
    integer(c_int),intent(in),optional :: idx(ndpth)
    integer*4 :: il,ierr
    integer(c_int) :: ndpthc,ltmtx_c,kt_c,nstrtl_c,ndl_c,nstrtt_c,ndt_c,nlf_c
    integer(c_int) :: idxc(ndpth+1)

    if(nlf==0) then
      st_lf=>null()
      return
    endif
    allocate(st_lf(nlf))
    if(ndpth>0.and.present(idx)) idxc(1:ndpth)=idx(1:ndpth)
    ndpthc=ndpth+1
    do il=1,nlf
      idxc(ndpthc)=il
      ierr=cHACApK_copy_leafmtx_c2f_scalars(ndpthc,idxc,ltmtx_c,kt_c,nstrtl_c,ndl_c,nstrtt_c,ndt_c,nlf_c)
      if(ierr/=0) return
      st_lf(il)%ltmtx=ltmtx_c
      st_lf(il)%kt=kt_c
      st_lf(il)%nstrtl=nstrtl_c
      st_lf(il)%ndl=ndl_c
      st_lf(il)%nstrtt=nstrtt_c
      st_lf(il)%ndt=ndt_c
      st_lf(il)%nlf=nlf_c
      if(nlf_c>0) then
        call HACApK_copy_leafmtx_c2f(st_lf(il)%st_lf,nlf_c,ndpthc,idxc)
      endif
    enddo
  end subroutine HACApK_copy_leafmtx_c2f

end module m_HACApK_base_c
