!=====================================================================*
!                                                                     *
!   Software Name : ppohBEM                                           *
!         Version : 0.100.0                                           *
!                                                                     *
!   License                                                           *
!     This file is part of ppohBEM.                                   *
!     ppohBEM is a free software, you can use it under the terms      *
!     of The MIT License (MIT). See LICENSE file and User's guide     *
!     for more details.                                               *
!                                                                     *
!   ppOpen-HPC project:                                               *
!     Open Source Infrastructure for Development and Execution of     *
!     Large-Scale Scientific Applications on Post-Peta-Scale          *
!     Supercomputers with Automatic Tuning (AT).                      *
!                                                                     *
!   Organizations:                                                    *
!     The University of Tokyo                                         *
!       - Information Technology Center                               *
!       - Atmosphere and Ocean Research Institute (AORI)              *
!       - Interfaculty Initiative in Information Studies              *
!         /Earthquake Research Institute (ERI)                        *
!       - Graduate School of Frontier Science                         *
!     Kyoto University                                                *
!       - Academic Center for Computing and Media Studies             *
!     Hokkaido University                                             *
!       - Information Initiative Center                               *
!     Japan Agency for Marine-Earth Science and Technology (JAMSTEC)  *
!                                                                     *
!   Sponsorship:                                                      *
!     Japan Science and Technology Agency (JST), Basic Research       *
!     Programs: CREST, Development of System Software Technologies    *
!     for post-Peta Scale High Performance Computing.                 *
!                                                                     *
!   Copyright (c) 2014 <Takeshi Iwashita, Takeshi Mifune, Yuki Noseda,*
!                    Yasuhito Takahashi, Masatoshi Kawai, Akihiro Ida>*
!                                                                     *
!=====================================================================*
!C***********************************************************************
!C  This file main program of Bem-BB,
!C  corrected the 'IF' block for PPOHBEM_INT_PARA_FC on Jan. 2017
!C  adjusted for Block Low Rank approximation on June 2017
!C  corrected format for VRK file on June 2017
!C  last modified by Akihiro Ida on June 2017,
!C***********************************************************************
!***ppohBEM_bem_bb_dense_mpi
program ppohBEM_bem_bb_dense_mpi
  use m_ppohBEM_bembb2hacapk
  implicit none
  include 'mpif.h'
  integer   irank, nrank
  integer   comm, ierr
  integer   iunit

  real*8  st1, st2, st21, st3, st4, time1, time2, time3, time4

  integer   ppohBEM_max_steps
  integer   ppohBEM_number_element_dof
  character ppohBEM_linear_solver*16
  real*8    ppohBEM_tor

  character filename*256

  integer   ppohBEM_nond, ppohBEM_nofc, ppohBEM_nond_on_face
  integer   ppohBEM_nint_para_fc, ppohBEM_ndble_para_fc

  type (coordinate), dimension(:), allocatable :: st_ppohBEM_np 

  integer, dimension(:,:), allocatable :: ppohBEM_face2node
  integer, dimension(:,:), allocatable :: ppohBEM_int_para_fc
  real*8,  dimension(:,:), allocatable :: ppohBEM_dble_para_fc

  integer  ndim, steps
  integer, EXTERNAL :: ppohBEM_pbicgstab_dense
  real*8, dimension(:,:), allocatable :: ppohBEM_a
  real*8, dimension(:),   allocatable :: ppohBEM_rhs
  real*8, dimension(:),   allocatable :: ppohBEM_sol

  integer  ext_ndim
  integer  lhp, ltp, proc_dim
  integer  act_nrank

  integer  i
  integer, dimension(:), allocatable :: act_nranks
  integer  MPI_GROUP_WORLD, act_grp
  integer  act_comm
  
!For HACApK
  character*32 value
  character*32 logfile
  real*8 ztol

  integer number,nlength,nstatus,nargs,il,it,nd,lrtrn
  type(st_HACApK_calc_entry) :: st_bemv
  type(st_HACApK_lcontrol) :: st_ctl
  logical ldense

  ldense=.false.
!  ldense=.true.

  comm = MPI_COMM_WORLD; ierr = 0

  call Initialization_MPI( comm, nrank, irank, ierr )

  if( ierr .ne. 0 ) then
    print*, 'Error: MPI_Init failed'
    goto 1000
  endif


  call Read_bem_bb_config( ppohBEM_number_element_dof, ppohBEM_linear_solver, &
                           ppohBEM_tor, ppohBEM_max_steps, irank, comm, ierr )

  if( ierr .ne. 0 ) then
    print*, 'Error: in reading bem_bb_config.txt'
    goto 1000
  endif

  call MPI_Barrier( comm, ierr )
  if( ierr .ne. 0 ) then
    print*, 'Error: MPI_Barrier #1'
    goto 1000
  endif

  st1 = MPI_Wtime()

  nargs=command_argument_count()
  if ( irank .eq. 0 ) then
    if(nargs==0)then
     filename =  "input_sample.pbf"; print*, 'Sample model data file is used.'
    else
      number=1
      call get_command_argument(number,value,nlength,nstatus)
      filename =trim(value)
    endif
  endif

  iunit = 10

  call Read_model_data( filename, ppohBEM_number_element_dof, &
                        ppohBEM_nond, ppohBEM_nofc, ppohBEM_nond_on_face, &
                        ppohBEM_nint_para_fc, ppohBEM_ndble_para_fc, &
                        st_ppohBEM_np, ppohBEM_face2node, &
                        ppohBEM_int_para_fc, ppohBEM_dble_para_fc, iunit, &
                        irank, comm, ierr )

  call MPI_Barrier( comm, ierr )
  if(ierr.ne.0) then
    print*, 'Error: MPI_Barrier #2'
    goto 1000
  endif

  st2 = MPI_Wtime()


  call Make_equation_data( ppohBEM_number_element_dof, &
                           ppohBEM_nond, ppohBEM_nofc, ppohBEM_nond_on_face, &
                           ppohBEM_nint_para_fc, ppohBEM_ndble_para_fc, &
                           st_ppohBEM_np, ppohBEM_face2node, &
                           ppohBEM_int_para_fc, ppohBEM_dble_para_fc, &
                           ndim, ppohBEM_a, ppohBEM_rhs, ppohBEM_sol, &
                           ext_ndim, lhp, ltp, proc_dim, irank, nrank, &
                           act_nrank,ldense, ierr )


  call MPI_Barrier( comm, ierr )
  if(ierr.ne.0) then
    print*, 'Error: MPI_Barrier #3'
    goto 1000
  endif

  st3 = MPI_Wtime()

    nd=ppohBEM_nofc*ppohBEM_number_element_dof
    lrtrn=Hacapk_Init(nd,st_ctl,st_bemv,comm)
    if(nargs==0)then
    else
      number=3
      call get_command_argument(number,value,nlength,nstatus)
      if(value=="") then; else
        read(value,*) st_ctl%param(42) ! BLR : leaf size for BLR 0:floor(sqrt(N))
      endif
      number=2
      call get_command_argument(number,value,nlength,nstatus)
      if(value=="") then; else
        read(value,*) st_ctl%param(21) ! cluster : leaf size
      endif
    endif
    
!    st_ctl%param(1)=0
    st_ctl%param(83)=ppohBEM_max_steps
    if    ( ppohBEM_linear_solver == "BICGSTAB" ) then
      st_ctl%param(85)=1
    elseif( ppohBEM_linear_solver == "GCRM" ) then
      st_ctl%param(85)=2
    else
      print*, 'Error: Invalid ppohBEM_linear_solver'
      goto 1000
    endif
!    if(st_ctl%param(1)>1)  print*,'func Hacapk_Init end'
    
    st_ctl%param(1) =2;        ! Print : 0:Only Error 1:STD 2:Dubug
    st_ctl%param(43)=1;       ! BH : Number of blocks with a color in a row of process grid
    st_ctl%param(53)=5;        ! H-matrix : maximun depth of block cluster tree 100
    st_ctl%param(60)=3;        ! 2:ACA+ 3:SVD 5:RRQR
    st_bemv%lperp=1            ! 0:scalar kernel, 1:vecter kernel
    
    lrtrn= bembb2hacapk(st_bemv, st_ctl, st_ppohBEM_np, ppohBEM_face2node,ppohBEM_dble_para_fc,ppohBEM_rhs,ppohBEM_sol,ppohBEM_tor, &
                        ppohBEM_int_para_fc,ppohBEM_nond,ppohBEM_nofc,ppohBEM_nond_on_face,ppohBEM_number_element_dof,&
                        ppohBEM_ndble_para_fc,ppohBEM_nint_para_fc)
  lrtrn=HACApK_finalize(st_ctl)
  
  call MPI_Barrier( comm, ierr )
  if(ierr.ne.0) then
    print*, 'Error: MPI_Barrier #22 !!!'
    goto 1000
  endif
  
  st21 = MPI_Wtime()
  
!  if(irank==0) print*,'H-matrix time=',st21-st3

if(ldense)then
  allocate( act_nranks(act_nrank), stat=ierr )
  if( ierr .ne. 0 ) then
    print*, 'Memory allocation #01 failed'
    goto 1000
  endif

  do i = 1, act_nrank
    act_nranks(i) = i - 1
  enddo

  call MPI_Comm_group( MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr )
  call MPI_Group_incl( MPI_GROUP_WORLD, act_nrank, act_nranks, act_grp, ierr )
  call MPI_Comm_create( MPI_COMM_WORLD, act_grp, act_comm, ierr )

  if( irank .lt. act_nrank ) then
    if( ppohBEM_linear_solver == "BICGSTAB" ) then
      steps = ppohBEM_pbicgstab_dense( irank, nrank, lhp, ltp, proc_dim, ndim,&
                                       ext_ndim, ppohBEM_a, ppohBEM_rhs, &
                                       ppohBEM_sol, ppohBEM_tor, &
                                       ppohBEM_max_steps, act_comm, ierr )
    endif
  endif

  call MPI_Comm_free( act_comm, ierr )
  call MPI_Group_free( MPI_GROUP_WORLD, ierr )
  call MPI_Group_free( act_grp, ierr )


  call MPI_Barrier( comm, ierr )
  if(ierr.ne.0) then
    print*, 'Error: MPI_Barrier #4'
    goto 1000
  endif

endif

  st4 = MPI_Wtime()

 if(.true.)then
  filename="sample.output.vtk"
  iunit = 11
  call Print_result( filename, ppohBEM_nond, ppohBEM_nofc, &
                     ppohBEM_nond_on_face, st_ppohBEM_np, ppohBEM_face2node, &
                     ndim, ppohBEM_sol, irank, iunit, ierr )
 endif

  time1 = st2-st1
  time3 = st21-st3
  time4 = st4-st1

  if ( irank .eq. 0 ) then
    print*,'before HACApK time = ',sngl(time1),' [sec]'
    print*,'       HACApK time = ',sngl(time3),' [sec]'
    print*,'_____________all time = ',sngl(time4),' [sec]'
  endif



1000  continue

  if(allocated(st_ppohBEM_np))         deallocate(st_ppohBEM_np)
  if(allocated(ppohBEM_face2node))     deallocate(ppohBEM_face2node)
  if(allocated(ppohBEM_int_para_fc))   deallocate(ppohBEM_int_para_fc)
  if(allocated(ppohBEM_dble_para_fc))  deallocate(ppohBEM_dble_para_fc)

  if(allocated(ppohBEM_a))             deallocate(ppohBEM_a)
  if(allocated(ppohBEM_rhs))           deallocate(ppohBEM_rhs)
  if(allocated(ppohBEM_sol))           deallocate(ppohBEM_sol)

  if(allocated(act_nranks))    deallocate(act_nranks)

  call MPI_Finalize (ierr)  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    subroutines    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!    Initialization_MPI    !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***Initialization_MPI
  subroutine Initialization_MPI( comm, nrank, irank, ierr )
    integer, intent(out) :: irank
    integer, intent(out) :: nrank
    integer, intent(in)  :: comm
    integer, intent(out) :: ierr

    call MPI_Init ( ierr )
    if( ierr .ne. 0 ) then
      print*, 'Error: MPI_Init failed !!!'
    endif

    call MPI_Comm_size ( comm, nrank, ierr )
    if( ierr .ne. 0 ) then
      print*, 'Error: MPI_Comm_size failed !!!'
    endif

    call MPI_Comm_rank ( comm, irank, ierr )
    if( ierr .ne. 0 ) then
      print*, 'Error: MPI_Comm_rank failed !!!'
    endif

    if (irank.eq.0) then
      print*, 'Number of processes ', nrank
    endif

    call MPI_Barrier( comm, ierr )
    if( ierr .ne. 0 ) then
      print*, 'Error: MPI_Barrier #1 !!!'
    endif

  end subroutine   !!!! subroutine Initialize_MPI


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!    Read_bem_bb_config    !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***Read_bem_bb_config
  subroutine Read_bem_bb_config( ppohBEM_number_element_dof, &
                                 ppohBEM_linear_solver, ppohBEM_tor, &
                                 ppohBEM_max_steps, irank, comm, ierr )

    integer,   intent(out) :: ppohBEM_number_element_dof
    character, intent(out) :: ppohBEM_linear_solver*16
    real*8,    intent(out) :: ppohBEM_tor
    integer,   intent(out) :: ppohBEM_max_steps
    integer,   intent(in)  :: irank
    integer,   intent(in)  :: comm
    integer,   intent(out) :: ierr

    if( irank .eq. 0 ) then
      open( 11, file="bem-bb-config.txt", action='read', iostat=ierr )
      read( 11, * ) ppohBEM_number_element_dof
      read( 11, * ) ppohBEM_linear_solver
      read( 11, * ) ppohBEM_tor
      read( 11, * ) ppohBEM_max_steps
      close( 11 )

      write(*,*) "Number of unknowns set on each face element = ", &
                  ppohBEM_number_element_dof
      write(*,*) "Selected linear solver = ", ppohBEM_linear_solver
      write(*,*) "Convergence criterion = ", ppohBEM_tor
      write(*,*) "Upper limit of iteration counts = ", ppohBEM_max_steps

    endif

    call MPI_Bcast ( ppohBEM_number_element_dof, 1, MPI_INTEGER, 0, comm, ierr )
    call MPI_Bcast ( ppohBEM_linear_solver, 16, MPI_CHARACTER, 0, comm, ierr )
    call MPI_Bcast ( ppohBEM_tor, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr )
    call MPI_Bcast ( ppohBEM_max_steps, 1, MPI_INTEGER, 0, comm, ierr )

  end subroutine   !!!! subroutine Read_bem_bb_config


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!    Read_model_data    !!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***Read_model_data
  subroutine Read_model_data( filename, ppohBEM_number_element_dof, &
                              ppohBEM_nond, ppohBEM_nofc, &
                              ppohBEM_nond_on_face, ppohBEM_nint_para_fc, &
                              ppohBEM_ndble_para_fc, &
                              st_ppohBEM_np, ppohBEM_face2node, &
                              ppohBEM_int_para_fc, ppohBEM_dble_para_fc, &
                              iunit, irank, comm, ierr )

    character, intent(in) :: filename*256
    integer,   intent(in)  :: ppohBEM_number_element_dof

    integer, intent(out) :: ppohBEM_nond
    integer, intent(out) :: ppohBEM_nofc
    integer, intent(out) :: ppohBEM_nond_on_face
    integer, intent(out) :: ppohBEM_nint_para_fc
    integer, intent(out) :: ppohBEM_ndble_para_fc

    type (coordinate), dimension(:), allocatable, intent(out) :: &
                                                              st_ppohBEM_np
    integer,           dimension(:,:), allocatable, intent(out) :: &
                                                              ppohBEM_face2node
    integer,           dimension(:,:), allocatable, intent(out) :: &
                                                            ppohBEM_int_para_fc
    real*8,            dimension(:,:), allocatable, intent(out) :: &
                                                           ppohBEM_dble_para_fc

    integer, intent(in)  :: iunit
    integer, intent(in)  :: irank
    integer, intent(in)  :: comm
    integer, intent(out) :: ierr

    integer  dims(5)

!    print*,'sub Read_model_data start'


    if ( irank .eq. 0 ) then
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Case: irank is equal to 0 (Master processes) !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call read_pbf( filename, &
                     ppohBEM_nond, ppohBEM_nofc, &
                     ppohBEM_nond_on_face, ppohBEM_nint_para_fc, &
                     ppohBEM_ndble_para_fc, &
                     st_ppohBEM_np, ppohBEM_face2node, &
                     ppohBEM_int_para_fc, ppohBEM_dble_para_fc, &
                     iunit, ierr )
      if( ierr .ne. 0 ) then
        call read_stl_ascii( filename, &
                             ppohBEM_nond, ppohBEM_nofc, &
                             ppohBEM_nond_on_face, ppohBEM_nint_para_fc, &
                             ppohBEM_ndble_para_fc, &
                             st_ppohBEM_np, ppohBEM_face2node, &
                             ppohBEM_int_para_fc, ppohBEM_dble_para_fc, &
                             iunit, ierr )
        if( ierr .ne. 0 ) then
          print*, 'Reading model data failed; exiting'
          stop
        endif
      endif

      dims(1) = ppohBEM_nond
      dims(2) = ppohBEM_nofc
      dims(3) = ppohBEM_nond_on_face
      dims(4) = ppohBEM_nint_para_fc
      dims(5) = ppohBEM_ndble_para_fc

      call MPI_Bcast( dims(1), 5, MPI_INTEGER, 0, comm, ierr )
      if( ierr .ne. 0 ) then
        print*, 'Error: MPI_Bcast#11'
      endif

      call MPI_Bcast( st_ppohBEM_np(1)%x, ppohBEM_nond*3, & 
                      MPI_DOUBLE_PRECISION, 0, comm, ierr )
      if( ierr .ne. 0 ) then
        print*, 'Error: MPI_Bcast#12'
      endif

      call MPI_Bcast ( ppohBEM_face2node(1,1), &
                       ppohBEM_nofc*ppohBEM_nond_on_face, MPI_INTEGER, 0, &
                       comm, ierr )
      if( ierr .ne. 0 ) then
        print*, 'Error: MPI_Bcast#13'
      endif

      if( ppohBEM_nint_para_fc .gt. 0 ) then
        call MPI_Bcast ( ppohBEM_int_para_fc(1,1), &
                         ppohBEM_nofc*ppohBEM_nint_para_fc, MPI_INTEGER, 0, &
                         comm, ierr )
        if( ierr .ne. 0 ) then
          print*, 'Error: MPI_Bcast#14 !!!'
        endif
      endif

      if( ppohBEM_ndble_para_fc .gt. 0 ) then
        call MPI_Bcast ( ppohBEM_dble_para_fc(1,1), &
                         ppohBEM_nofc*ppohBEM_ndble_para_fc, &
                         MPI_DOUBLE_PRECISION, 0, comm, ierr )
        if( ierr .ne. 0 ) then
          print*, 'Error: MPI_Bcast#15 !!!'
        endif
      endif

      close ( iunit )

    else
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Case: irank is not equal to 0 (Worker processes) !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call MPI_Bcast( dims(1), 5, MPI_INTEGER, 0, comm, ierr )
      if( ierr .ne. 0 ) then
        print*, 'Error: MPI_Bcast#21 !!!'
      endif

      ppohBEM_nond = dims(1)
      ppohBEM_nofc = dims(2)
      ppohBEM_nond_on_face  = dims(3)
      ppohBEM_nint_para_fc  = dims(4)
      ppohBEM_ndble_para_fc = dims(5)

      allocate( st_ppohBEM_np(ppohBEM_nond), stat=ierr )
      if( ierr .ne. 0 ) then
        print*, 'Memory allocation #21 failed !!!'
      endif
      call MPI_Bcast( st_ppohBEM_np(1)%x, ppohBEM_nond*3, &
                      MPI_DOUBLE_PRECISION, 0, comm, ierr )
      if( ierr .ne. 0 ) then
        print*, 'Error: MPI_Bcast#22'
      endif

      allocate( ppohBEM_face2node(ppohBEM_nond_on_face, ppohBEM_nofc), &
                stat=ierr )
      if( ierr .ne. 0 ) then
        print*, 'Memory allocation #22 failed !!!'
      endif
      call MPI_Bcast ( ppohBEM_face2node(1,1), &
                       ppohBEM_nofc*ppohBEM_nond_on_face, MPI_INTEGER, 0, &
                       comm, ierr )
      if( ierr .ne. 0 ) then
        print*, 'Error: MPI_Bcast#23'
      endif

      if ( ppohBEM_nint_para_fc .gt. 0 ) then
        allocate( ppohBEM_int_para_fc(ppohBEM_nint_para_fc, ppohBEM_nofc), &
                  stat=ierr )
        if( ierr .ne. 0 ) then
          print*, 'Memory allocation #23 failed !!!'
        endif
        call MPI_Bcast ( ppohBEM_int_para_fc(1,1), &
                         ppohBEM_nofc*ppohBEM_nint_para_fc, MPI_INTEGER, 0, &
                         comm, ierr )
        if( ierr .ne. 0 ) then
          print*, 'Error: MPI_Bcast#24 !!!'
        endif
      else
        allocate( ppohBEM_int_para_fc(1,1), stat=ierr )
      endif
      if ( ppohBEM_ndble_para_fc .gt. 0 ) then
        allocate( ppohBEM_dble_para_fc(ppohBEM_ndble_para_fc, ppohBEM_nofc), & 
                  stat=ierr )
        if( ierr .ne. 0 ) then
          print*, 'Memory allocation #24 failed !!!'
        endif
        call MPI_Bcast ( ppohBEM_dble_para_fc(1,1), &
                         ppohBEM_nofc*ppohBEM_ndble_para_fc, &
                         MPI_DOUBLE_PRECISION, 0, comm, ierr )
        if( ierr .ne. 0 ) then
          print*, 'Error: MPI_Bcast#25 !!!'
        endif
      else
        allocate( ppohBEM_dble_para_fc(1, 1),stat=ierr )
      endif
    endif
!    print*,'sub Read_model_data end'
    
  end subroutine   !!!! subroutine Read_model_data


  subroutine read_pbf( filename, &
                       ppohBEM_nond, ppohBEM_nofc, &
                       ppohBEM_nond_on_face, ppohBEM_nint_para_fc, &
                       ppohBEM_ndble_para_fc, &
                       st_ppohBEM_np, ppohBEM_face2node, &
                       ppohBEM_int_para_fc, ppohBEM_dble_para_fc, &
                       iunit, ierr )

    character, intent(in) :: filename*256

    integer, intent(out) :: ppohBEM_nond
    integer, intent(out) :: ppohBEM_nofc
    integer, intent(out) :: ppohBEM_nond_on_face
    integer, intent(out) :: ppohBEM_nint_para_fc
    integer, intent(out) :: ppohBEM_ndble_para_fc

    type (coordinate), dimension(:), allocatable, intent(out) :: &
                                                              st_ppohBEM_np
    integer,           dimension(:,:), allocatable, intent(out) :: &
                                                              ppohBEM_face2node
    integer,           dimension(:,:), allocatable, intent(out) :: &
                                                            ppohBEM_int_para_fc
    real*8,            dimension(:,:), allocatable, intent(out) :: &
                                                           ppohBEM_dble_para_fc

    integer, intent(in)  :: iunit
    integer, intent(out) :: ierr

    integer  i, j

    open( iunit, file=filename, action='read', pad='yes', iostat=ierr )
    if( ierr .ne. 0 ) then
      print*, trim(filename), ' does not exist'; stop
    endif

    !!!!  Read number of nodes from input data file : ppohBEM_nond  !!!!
    read( iunit, *, iostat=ierr ) ppohBEM_nond
    if( ierr .ne. 0 ) then
      print*, 'Input file is not in pbf format'
      return
    endif

    !!!!  Allocation for the array for the coordinates of the nodes : 
    !!!!      st_ppohBEM_np !!!!
    allocate( st_ppohBEM_np(ppohBEM_nond), stat=ierr )
    if( ierr .ne. 0 ) then
      print*, 'Memory allocation #1 failed'
    endif

    !!!!  Read the coordinates of the nodes from input data file :
    !!!!      st_ppohBEM_np  !!!!
    do i = 1, ppohBEM_nond
      read( iunit, * ) st_ppohBEM_np(i)%x, st_ppohBEM_np(i)%y, &
          st_ppohBEM_np(i)%z
    enddo

    !!!!  Read number of faces from input data file : ppohBEM_nofc  !!!!
    read( iunit, * ) ppohBEM_nofc

    !!!!  Read number of nodes on each face from input data file : 
    !!!!       ppohBEM_nond_on_face  !!!!

    read( iunit, * ) ppohBEM_nond_on_face

    !!!!  Read number of integer parameters set on each face from input data
    !!!!    file : ppohBEM_nint_para_fc  !!!!
          
    read( iunit, * ) ppohBEM_nint_para_fc

    !!!!  Read number of double precision parameters set on each face from
    !!!!    input data file : ppohBEM_ndble_para_fc  !!!!
    
    read( iunit, * ) ppohBEM_ndble_para_fc

    allocate( ppohBEM_face2node(ppohBEM_nond_on_face, ppohBEM_nofc), &
              stat=ierr )
    if( ierr .ne. 0 ) then
      print*, 'Memory allocation #11 failed'
    endif
    do i = 1, ppohBEM_nofc
      read( iunit, * ) (ppohBEM_face2node(j,i), j=1,ppohBEM_nond_on_face)
    enddo

    if( ppohBEM_nint_para_fc .gt. 0 ) then
      allocate( ppohBEM_int_para_fc(ppohBEM_nint_para_fc, ppohBEM_nofc), &
                stat=ierr )
      if( ierr .ne. 0 ) then
        print*, 'Memory allocation #12 failed'
      endif
      do i = 1, ppohBEM_nofc
        read(iunit, * ) (ppohBEM_int_para_fc(j,i), j=1, ppohBEM_nint_para_fc)
      enddo
    else
      allocate( ppohBEM_int_para_fc(1,1), stat=ierr )
    endif

    if( ppohBEM_ndble_para_fc .gt. 0 ) then
      allocate( ppohBEM_dble_para_fc(ppohBEM_ndble_para_fc, ppohBEM_nofc), &
                stat=ierr )
      if( ierr .ne. 0 ) then
        print*, 'Memory allocation #13 failed'
      endif
      do i = 1, ppohBEM_nofc
        read( iunit, * ) ( ppohBEM_dble_para_fc(j,i), &
                           j=1, ppohBEM_ndble_para_fc)
      enddo
    else
      allocate( ppohBEM_dble_para_fc(1, 1),stat=ierr )
    endif

  end subroutine   !!! subroutine read_pbf


  subroutine read_stl_ascii( filename, &
                             ppohBEM_nond, ppohBEM_nofc, &
                             ppohBEM_nond_on_face, ppohBEM_nint_para_fc, &
                             ppohBEM_ndble_para_fc, &
                             st_ppohBEM_np, ppohBEM_face2node, &
                             ppohBEM_int_para_fc, ppohBEM_dble_para_fc, &
                             iunit, ierr )

    use m_stl, only: t_stl, STL_Read_ascii, STL_Remove_duplicate_vertex, STL_Finalize

    character, intent(in) :: filename*256

    integer, intent(out) :: ppohBEM_nond
    integer, intent(out) :: ppohBEM_nofc
    integer, intent(out) :: ppohBEM_nond_on_face
    integer, intent(out) :: ppohBEM_nint_para_fc
    integer, intent(out) :: ppohBEM_ndble_para_fc

    type (coordinate), dimension(:), allocatable, intent(out) :: &
                                                              st_ppohBEM_np
    integer,           dimension(:,:), allocatable, intent(out) :: &
                                                              ppohBEM_face2node
    integer,           dimension(:,:), allocatable, intent(out) :: &
                                                            ppohBEM_int_para_fc
    real*8,            dimension(:,:), allocatable, intent(out) :: &
                                                           ppohBEM_dble_para_fc

    integer, intent(in)  :: iunit
    integer, intent(out) :: ierr

    type(t_stl) :: stl
    integer  i, j

    call STL_Read_ascii( stl, filename, ierr )
    if( ierr .ne. 0 ) then
      print*, 'Input file is not in STL ascii format'
      return
    endif
    call STL_Remove_duplicate_vertex( stl )

    ppohBEM_nond = stl%nvertex
    ppohBEM_nofc = stl%nfacet
    ppohBEM_nond_on_face = 3
    ppohBEM_nint_para_fc = 0
    ppohBEM_ndble_para_fc = 0

    allocate( st_ppohBEM_np(ppohBEM_nond), stat=ierr )
    if( ierr .ne. 0 ) then
      print*, 'Memory allocation #1 failed'
    endif

    do i = 1, ppohBEM_nond
      st_ppohBEM_np(i)%x = stl%vertex(1,i)
      st_ppohBEM_np(i)%y = stl%vertex(2,i)
      st_ppohBEM_np(i)%z = stl%vertex(3,i)
    enddo

    allocate( ppohBEM_face2node(ppohBEM_nond_on_face, ppohBEM_nofc), &
              stat=ierr )
    if( ierr .ne. 0 ) then
      print*, 'Memory allocation #11 failed'
    endif
    do i = 1, ppohBEM_nofc
      do j = 1, 3
        ppohBEM_face2node(j,i) = stl%v_on_f(j,i) - 1
      enddo
    enddo

    allocate( ppohBEM_int_para_fc(1,1), stat=ierr )

    allocate( ppohBEM_dble_para_fc(1, 1),stat=ierr )

    call STL_Finalize( stl )

  end subroutine   !!! subroutine read_stl_ascii


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!    Make_equation_data    !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***Make_equation_data
  subroutine Make_equation_data( ppohBEM_number_element_dof, &
                                 ppohBEM_nond, ppohBEM_nofc, &
                                 ppohBEM_nond_on_face, ppohBEM_nint_para_fc, &
                                 ppohBEM_ndble_para_fc, &
                                 st_ppohBEM_np, ppohBEM_face2node, &
                                 ppohBEM_int_para_fc, ppohBEM_dble_para_fc, &
                                 ndim, ppohBEM_a, ppohBEM_rhs, ppohBEM_sol, &
                                 ext_ndim, lhp, ltp, proc_dim, irank, nrank, &
                                 act_nrank,ldense, ierr )
  use m_ppohBEM_bembb2hacapk

    integer, intent(in)  :: ppohBEM_number_element_dof

    integer, intent(in)  :: ppohBEM_nond
    integer, intent(in)  :: ppohBEM_nofc
    integer, intent(in)  :: ppohBEM_nond_on_face
    integer, intent(in)  :: ppohBEM_nint_para_fc
    integer, intent(in)  :: ppohBEM_ndble_para_fc

    type (coordinate), dimension(:),  allocatable, intent(in) :: st_ppohBEM_np
    integer, dimension(:,:), allocatable, intent(in) :: ppohBEM_face2node
    integer, dimension(:,:), allocatable, intent(in) :: ppohBEM_int_para_fc
    real*8, dimension(:,:), allocatable, intent(in) :: ppohBEM_dble_para_fc

    integer, intent(out) :: ndim
    real*8, dimension(:,:), allocatable, intent(out) :: ppohBEM_a
    real*8, dimension(:),   allocatable, intent(out) :: ppohBEM_rhs
    real*8, dimension(:),   allocatable, intent(out) :: ppohBEM_sol

    integer, intent(out) :: ext_ndim
    integer, intent(out) :: lhp, ltp
    integer, intent(out) :: proc_dim
    integer, intent(in)  :: irank
    integer, intent(in)  :: nrank
    integer, intent(out) :: act_nrank
    integer, intent(out) :: ierr

    integer  i, j

    integer  thr_dim, thr_dim1
    integer  ith, nth, mth
    integer  j_st, j_en

    logical ldense

!    print*,'sub Make_equation_data start'

    ndim = ppohBEM_nofc * ppohBEM_number_element_dof

    proc_dim = ndim / nrank

    if( ndim .eq. proc_dim * nrank ) then
      lhp = irank * proc_dim
      act_nrank = nrank
    else
      proc_dim = proc_dim + 1
      lhp = irank * proc_dim
      act_nrank = (ndim / proc_dim) + 1
    endif

    ext_ndim = act_nrank * proc_dim

    ltp = lhp + proc_dim
    lhp = lhp + 1

    if(ldense) allocate( ppohBEM_a(ext_ndim,lhp:ltp), stat = ierr )
    allocate( ppohBEM_rhs(ext_ndim), stat = ierr )
    allocate( ppohBEM_sol(ext_ndim), stat = ierr )


!!!!!!!!!!!! Thread parallelization starts !!!!!!!!!!!!
!$omp parallel private(thr_dim, thr_dim1, ith, nth, mth, j_st, j_en, i, j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ith = omp_get_thread_num()
    nth = omp_get_num_threads()

    if ( ith.eq.0 .and. irank.eq.0 ) then
      print*, 'Number of total threads, Thread number ', nth, ith
      print*, 'ext_ndim = ', ext_ndim, ', ndim = ', ndim
    endif

    thr_dim = ndim / nth
    mth  = ndim - thr_dim * nth
    thr_dim1 = thr_dim + 1
    j_st = 0
    if( ith .lt. mth ) then
      do i = 1, ith
        j_st = j_st + thr_dim1
      enddo
      j_en = j_st + thr_dim1
    else
      do i = 1, mth
        j_st = j_st + thr_dim1
      enddo
      do i = mth+1, ith
        j_st = j_st + thr_dim
      enddo
      j_en = j_st + thr_dim
    endif

  if(ldense)then
  !$omp do private(j) schedule(static)
    do i = lhp, ltp
      do j = 1, ext_ndim
        ppohBEM_a(j,i) = 0.0d0
      enddo
    enddo
  !$omp end do
  endif
    if( ltp .gt. ndim ) then
      ltp = ndim
    endif

    j_st = j_st + 1


!!!!!!!!!!!!!!! Caution: the function element_ij should be given by users. !!!!
  if(ldense)then
    do i=lhp, ltp
      do j=j_st, j_en
        ppohBEM_a(j,i) = ppohBEM_matrix_element_ij( i, j, ppohBEM_nond, &
                                     ppohBEM_nofc, &
                                     ppohBEM_nond_on_face, st_ppohBEM_np, &
                                     ppohBEM_int_para_fc, &
                                     ppohBEM_nint_para_fc, &
                                     ppohBEM_dble_para_fc, &
                                     ppohBEM_ndble_para_fc, &
                                     ppohBEM_face2node,st_bemv%lperp )
      enddo
    enddo
  endif
    do j=j_st, j_en
      ppohBEM_rhs(j) = ppohBEM_right_hand_side_vector_element_i(j, &
                            ppohBEM_nond, &
                            ppohBEM_nofc, ppohBEM_nond_on_face, &
                            st_ppohBEM_np, ppohBEM_int_para_fc, &
                            ppohBEM_nint_para_fc, ppohBEM_dble_para_fc, &
                            ppohBEM_ndble_para_fc, ppohBEM_face2node )

      ppohBEM_sol(j) = 0.0d0
    enddo

!!!!!!!!!!!!!!!!!!
!$omp end parallel
!!!!!!!!!!!!!!!!!!
!    print*,'sub Make_equation_data end'

  end subroutine   !!!! subroutine Make_equation_data


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!    Print_result    !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***Print_result
  subroutine Print_result( filename, ppohBEM_nond, ppohBEM_nofc, &
                           ppohBEM_nond_on_face, st_ppohBEM_np, &
                           ppohBEM_face2node, ndim, ppohBEM_sol, &
                           irank, iunit, ierr )

    character, intent(in) :: filename*256
    integer,   intent(in)  :: ppohBEM_nond
    integer,   intent(in)  :: ppohBEM_nofc
    integer,   intent(in)  :: ppohBEM_nond_on_face

    type (coordinate), dimension(:),   allocatable, intent(in) :: st_ppohBEM_np
    integer, dimension(:,:), allocatable, intent(in) :: ppohBEM_face2node

    integer, intent(in) :: ndim
    real*8, dimension(:),   allocatable, intent(in) :: ppohBEM_sol

    integer,  intent(in)  :: irank
    integer,  intent(in)  :: iunit
    integer,  intent(out) :: ierr

    integer  i, j

    ierr=0

    if( irank .ne. 0 ) then
      goto 5000
    endif

    open( iunit, file=filename, status='replace')
    write( iunit, '(a)' ) "# vtk DataFile Version 2.0"
    write( iunit, '(a)' ) "header : this data is results of BEM-BB"
    write( iunit, '(a)' ) "ASCII"
    write( iunit, '(a)' ) "DATASET UNSTRUCTURED_GRID"
    write( iunit, '(a,i0,a)' ) "POINTS ", ppohBEM_nond, " float"
    do i = 1, ppohBEM_nond
!      write( iunit, '(f16.12, 1x, f16.12, 1x, f16.12)' ) &
      write( iunit, '(1pe15.7, 1x, 1pe15.7, 1x, 1pe15.7)' ) &
             st_ppohBEM_np(i)%x, st_ppohBEM_np(i)%y, st_ppohBEM_np(i)%z
    enddo
    write( iunit, '(a,i0,1x,i0)' ) "CELLS ", ppohBEM_nofc, &
                                   (ppohBEM_nond_on_face+1)*ppohBEM_nofc
    do i = 1, ppohBEM_nofc
      write( iunit, '(i0,1x)', advance='no' ) ppohBEM_nond_on_face
      do j = 1, ppohBEM_nond_on_face
        write( iunit, '(i0,1x)', advance='no' ) ppohBEM_face2node(j,i)
      enddo
      write( iunit, * ) " "
    enddo
    write( iunit, '(a,i0)' ) "CELL_TYPES ", ppohBEM_nofc
    do i = 1, ppohBEM_nofc
      write( iunit, * ) "5"
    enddo
    write( iunit, * ) "CELL_DATA ", ppohBEM_nofc
    write( iunit, * ) "SCALARS solve float 1"
    write( iunit, * ) "LOOKUP_TABLE default"
    do i = 1, ppohBEM_nofc
      write( iunit, '(f32.28)' ) ppohBEM_sol(i)
    enddo

    close( iunit )

5000 continue
    return
  end subroutine   !!!! subroutine Print_result

end  !!!!!!!!!!!!! program ppohBEM_bem_bb_dense_mpi !!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!    pbicgstab    !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***ppohBEM_pbicgstab_dense
integer function ppohBEM_pbicgstab_dense( irank, nrank, lhp, ltp, proc_dim, &
                                          ndim, ext_ndim, a, rhs, sol, &
                                          ppohBEM_tor, ppohBEM_max_steps, &
                                          act_comm, ierr )
  implicit none
  integer omp_get_thread_num, omp_get_num_threads

  include 'mpif.h'

  integer, intent(in) :: irank, nrank
  integer, intent(in) :: lhp, ltp
  integer, intent(in) :: proc_dim, ndim, ext_ndim

  real*8, dimension( ext_ndim, lhp:ltp ), intent(in) :: a
  real*8, dimension( ext_ndim ),          intent(in) :: rhs
  real*8, dimension( ext_ndim ),          intent(inout) :: sol

  real*8,  intent(in) :: ppohBEM_tor
  integer, intent(in) :: ppohBEM_max_steps

  integer, intent(in) :: act_comm
  integer, intent(out) :: ierr

  real*8, dimension(:), allocatable :: r, shdw, p, t, kp, akp, kt, akt
  real*8, dimension(:), allocatable :: asum
  real*8   asum1, asum2

  integer  step, nsteps, i, j
  real*8   alpha, beta, zeta, nom, den, nomold, rnorm, bnorm

  integer  thr_dim, thr_dim1
  integer  ith, nth, inth, mth
  integer  j_st, j_en, i_st, i_en

  integer  iunit

!!!! Allocation of arrays !!!!
  allocate( r(ext_ndim), shdw(ext_ndim), p(ext_ndim), t(ext_ndim), &
            kp(ext_ndim), stat=ierr )
  if( ierr .ne. 0 ) then
    print*, 'Memory allocation #101 failed !!!'
    goto 1000
  endif
  allocate( akp(ext_ndim), kt(ext_ndim), akt(ext_ndim), stat=ierr)
  if( ierr .ne. 0 ) then
    print*, 'Memory allocation #102 failed !!!'
    goto 1000
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp parallel &
!$omp private(step,i,j,alpha,beta,zeta,nom,den,nomold,rnorm,bnorm,thr_dim,thr_dim1,ith,nth,inth,mth,j_st,j_en,i_st,i_en)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ith = omp_get_thread_num()
  nth = omp_get_num_threads()

  thr_dim  = ndim / nth
  mth = ndim - thr_dim * nth
  thr_dim1 = thr_dim + 1
  j_st = 0
  if( ith .lt. mth ) then
    do i = 1, ith
      j_st = j_st + thr_dim1
    enddo
    j_en = j_st + thr_dim1
  else
    do i = 1, mth
      j_st = j_st + thr_dim1
    enddo
    do i = mth+1, ith
      j_st = j_st + thr_dim
    enddo
    j_en = j_st + thr_dim
  endif

  thr_dim = proc_dim / nth
  mth  = proc_dim - thr_dim * nth
  thr_dim1 = thr_dim + 1
  i_st = lhp - 1
  if( ith .lt. mth ) then
    do i = 1, ith
      i_st = i_st + thr_dim1
    enddo
    i_en = i_st + thr_dim1
  else
    do i = 1, mth
      i_st = i_st + thr_dim1
    enddo
    do i = mth+1, ith
      i_st = i_st + thr_dim
    enddo
    i_en = i_st + thr_dim
  endif

  j_st = j_st + 1
  i_st = i_st + 1

  ith  = ith + 1
  inth = ith + nth


!!!! Initialization !!!!
  do i = i_st, i_en
    r(i) = 0.0d0;   shdw(i) = 0.0d0
    t(i) = 0.0d0;   akp(i) = 0.0d0
    kt(i) = 0.0d0;  akt(i) = 0.0d0
  enddo
  do j = j_st, j_en
    r(j) = 0.0d0;   shdw(j) = 0.0d0
    p(j) = 0.0d0;   t(j) = 0.0d0
    kp(j) = 0.0d0;  akp(j) = 0.0d0
    kt(j) = 0.0d0;  akt(j) = 0.0d0
  enddo
  if( j_en==ndim ) then
    do j = ndim, ext_ndim
      r(j) = 0.0d0;   shdw(j) = 0.0d0
      p(j) = 0.0d0;   t(j) = 0.0d0

      kp(j) = 0.0d0;  akp(j) = 0.0d0
      kt(j) = 0.0d0;  akt(j) = 0.0d0
    enddo
  endif
  alpha = 0.0d0;  beta = 0.0d0;  zeta = 0.0d0



!$omp master
  allocate( asum(nth+nth), stat=ierr )
  asum(1) = 0.0d0
!$omp end master
  if( ierr .ne. 0 ) then
    print*, 'Memory allocation #103 failed !!!'
    goto 990
  endif
!$omp barrier

  call residual_direct( i_st, i_en, ndim, ext_ndim, lhp, ltp, a, sol, rhs, r )

!$omp master
  call MPI_Allgather ( MPI_IN_PLACE, proc_dim, MPI_DOUBLE_PRECISION, &
                       r, proc_dim, MPI_DOUBLE_PRECISION, act_comm, ierr )
!$omp end master
!$omp barrier


!!!! Set shadow vector !!!!
  do j = j_st, j_en
    shdw(j) = r(j)
  enddo


!$omp barrier
  call dot_product( i_st, i_en, ext_ndim, rhs, rhs, asum(ith) )
  call dot_product( i_st, i_en, ext_ndim, r, r, asum(inth) )
  asum1 = 0.0
  asum2 = 0.0
!$omp barrier
!$omp do reduction(+:asum1, asum2) schedule(static,1)
  do i = 1, nth
    asum1 = asum1 + asum(i)
    asum2 = asum2 + asum(i+nth)
  enddo
!$omp barrier
!$omp master
  asum(1) = asum1
  asum(2) = asum2
  call MPI_Allreduce( MPI_IN_PLACE, asum, 2, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      act_comm, ierr )
!$omp end master
!$omp barrier
  bnorm = sqrt( asum(1) )
  rnorm = sqrt( asum(2) )
  asum1 = 0.0
  asum2 = 0.0
!$omp barrier


  if( irank .eq. 0 ) then
!$omp master
    write(*,1010) 'Original relative residual norm =', rnorm/bnorm
!$omp end master
  endif

  if ( rnorm .lt. ppohBEM_tor*bnorm ) then
    nsteps = 0; step = 0;  goto 900
  endif

!!!! BiCGSTAB iteration !!!!
  do step=1, ppohBEM_max_steps

    do j = j_st, j_en
      p(j) = r(j) + beta * ( p(j) - zeta * akp(j) )
    enddo
  !!!! No preconditioning !!!!
    do j = j_st, j_en
      kp(j) = p(j)
    enddo
  !$omp barrier
    call matvec_direct( i_st, i_en, ndim, ext_ndim, lhp, ltp, a, kp, akp )
  !$omp barrier
  !$omp master
    call MPI_Allgather( MPI_IN_PLACE, proc_dim, MPI_DOUBLE_PRECISION, &
                        akp, proc_dim, MPI_DOUBLE_PRECISION, act_comm, ierr )
  !$omp end master
  !$omp barrier
    call dot_product( i_st, i_en, ext_ndim, shdw, r, asum(ith) )
    call dot_product( i_st, i_en, ext_ndim, shdw, akp, asum(inth) )
  !$omp barrier
  !$omp do reduction(+:asum1, asum2) schedule(static,1)
    do i = 1, nth
      asum1 = asum1 + asum(i)
      asum2 = asum2 + asum(i+nth)
    enddo
  !$omp barrier
  !$omp master
    asum(1) = asum1
    asum(2) = asum2
    call MPI_Allreduce( MPI_IN_PLACE, asum, 2, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, act_comm, ierr )
  !$omp end master
  !$omp barrier
    nom = asum(1)
    den = asum(2)
    asum1 = 0.0
    asum2 = 0.0
  !$omp barrier
    alpha  = nom / den
    nomold = nom

    do j = j_st, j_en
      t(j) = r(j) - alpha * akp(j)
    enddo

    !!!! No preconditioning !!!!
    do j = j_st, j_en
      kt(j) = t(j)
    enddo

  !$omp barrier
    call matvec_direct( i_st, i_en, ndim, ext_ndim, lhp, ltp, a, kt, akt )
  !$omp barrier
  !$omp master
    call MPI_Allgather( MPI_IN_PLACE, proc_dim, MPI_DOUBLE_PRECISION, &
                        akt, proc_dim, MPI_DOUBLE_PRECISION, act_comm, ierr )
  !$omp end master
  !$omp barrier
    call dot_product( i_st, i_en, ext_ndim, akt, t, asum(ith) )
    call dot_product( i_st, i_en, ext_ndim, akt, akt, asum(inth) )
  !$omp barrier
  !$omp do reduction(+:asum1, asum2) schedule(static,1)
    do i = 1, nth
      asum1 = asum1 + asum(i)
      asum2 = asum2 + asum(i+nth)
    enddo
  !$omp barrier
  !$omp master
    asum(1) = asum1
    asum(2) = asum2
    call MPI_Allreduce( MPI_IN_PLACE, asum, 2, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        act_comm, ierr )
  !$omp end master
  !$omp barrier
    nom = asum(1)
    den = asum(2)
    asum1 = 0.0
    asum2 = 0.0
  !$omp barrier
    zeta  = nom / den

    do j = j_st, j_en
      sol(j) = sol(j) + alpha * kp(j) + zeta * kt(j)
    enddo
    do j = j_st, j_en
      r(j) = t(j) - zeta * akt(j)
    enddo

  !$omp barrier
    call dot_product( i_st, i_en, ext_ndim, shdw, r, asum(ith) )
    call dot_product( i_st, i_en, ext_ndim, r, r, asum(inth) )
  !$omp barrier
  !$omp do reduction(+:asum1, asum2) schedule(static,1)
    do i = 1, nth
      asum1 = asum1 + asum(i)
      asum2 = asum2 + asum(i+nth)
    enddo
  !$omp barrier
  !$omp master
    asum(1) = asum1
    asum(2) = asum2
    call MPI_Allreduce( MPI_IN_PLACE, asum, 2, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, act_comm, ierr )
  !$omp end master
  !$omp barrier
    beta = alpha * asum(1) / ( zeta * nomold )
    rnorm = sqrt( asum(2) )
    asum1 = 0.0
    asum2 = 0.0
  !$omp barrier
    if( irank .eq. 0 ) then
  !$omp master
      write(*,1011) 'Step', step, ' relative residual =', rnorm/bnorm
  !$omp end master
    endif
    if( rnorm .lt. ppohBEM_tor*bnorm ) then
      exit   !!! end iteration
    endif
  enddo  !!!! do step=1, ppohBEM_max_steps !!!!

900  continue

  call residual_direct( i_st, i_en, ndim, ext_ndim, lhp, ltp, a, sol, rhs, r )

!$omp master
  call MPI_Allgather ( MPI_IN_PLACE, proc_dim, MPI_DOUBLE_PRECISION, &
                       r, proc_dim, MPI_DOUBLE_PRECISION, act_comm, ierr )
!$omp end master
!$omp barrier

  call dot_product( i_st, i_en, ext_ndim, r, r, asum(ith) )

!$omp barrier
!$omp do reduction(+:asum1) schedule(static,1)
  do i = 1, nth
    asum1 = asum1 + asum(i)
  enddo
!$omp barrier
!$omp master
  asum(1) = asum1
  call MPI_Allreduce( MPI_IN_PLACE, asum, 1, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, act_comm, ierr )
!$omp end master
!$omp barrier
  rnorm = sqrt( asum(1) )
!$omp barrier

    if( irank .eq. 0 ) then
!$omp master
      write(*,*) ' '
      write(*,1010) 'Relative residual norm = ', rnorm/bnorm
      write(*,*) ' '
      write(*,*) ' '
!$omp end master
    endif
    nsteps = step

990  continue

!$omp end parallel

1000 continue
  if(allocated(r))    deallocate(r)
  if(allocated(shdw)) deallocate(shdw)
  if(allocated(p))    deallocate(p)
  if(allocated(t))    deallocate(t)
  if(allocated(kp))   deallocate(kp)
  if(allocated(akp))  deallocate(akp)
  if(allocated(kt))   deallocate(kt)
  if(allocated(akt))  deallocate(akt)
  if(allocated(asum)) deallocate(asum)

  ppohBEM_pbicgstab_dense = nsteps
  return

!1010 FORMAT(a, 1ES)
!1011 FORMAT(a, 1I5, a, 1ES)

1010 FORMAT(a, 1pe15.7)
1011 FORMAT(a, 1I5, a, 1pe15.7)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    subroutines    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!    residual_direct    !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***residual_direct
  subroutine residual_direct( i_st, i_en, ndim, ext_ndim, lhp, ltp, a, &
                              x, rhs, r )

    integer, intent(in) :: i_st, i_en, ndim
    integer, intent(in) :: ext_ndim, lhp, ltp

    real*8, dimension( ext_ndim, lhp:ltp ), intent(in) :: a
    real*8, dimension( ext_ndim ),          intent(in) :: x, rhs
    real*8, dimension( ext_ndim ),          intent(out) :: r

    integer  i, j

    do i = i_st, i_en
      r(i) = rhs(i)
    enddo

    do i = i_st, i_en
      do j = 1, ndim
        r(i) = r(i) - a(j,i) * x(j)
      enddo
    enddo

  end subroutine   !!!! subroutine residual_direct

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!    matvec_direct    !!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***matvec_direct
  subroutine matvec_direct( i_st, i_en, ndim, ext_ndim, lhp, ltp, a, p, q )

    integer, intent(in) :: i_st, i_en, ndim
    integer, intent(in) :: ext_ndim, lhp, ltp

    real*8, dimension( ext_ndim, lhp:ltp ), intent(in) :: a
    real*8, dimension( ext_ndim ),          intent(in) :: p
    real*8, dimension( ext_ndim ),          intent(out) :: q

    integer  i, j

    do i = i_st, i_en
      q(i) = 0.0d0
    enddo

    do i = i_st, i_en
      do j = 1, ndim
        q(i) = q(i) + a(j,i) * p(j)
      enddo
    enddo

  end subroutine   !!!! subroutine residual_direct

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!    dot_product    !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***dot_product
  subroutine dot_product( i_st, i_en, ext_ndim, p, q, asum )

    integer, intent(in) :: i_st, i_en, ext_ndim

    real*8, dimension( ext_ndim ), intent(in) :: p, q
    real*8, intent(out) :: asum

    integer  i

    asum = 0.0d0
    do i = i_st, i_en
      asum = asum + p(i) * q(i)
    enddo

  end subroutine   !!!! subroutine dot_product

end function ppohBEM_pbicgstab_dense

