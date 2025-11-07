!=====================================================================*
!                                                                     *
!   Software Name : ppohBEM                                           *
!         Version : 0.1                                               *
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
module m_ppohBEM_bembb2hacapk
 use m_HACApK_use
 implicit real*8(a-h,o-z)
contains

!*** bembb2hacapk
! np                ; coordinates of vertexes: x,y,z
! face2node         ; nodes which compose face
! dble_para_fc      ; real*8 parameters sat on each face
! rhs               ; right hand side vector
! sol               ; solution vector
! param             ; parameter for hacapk
! lpmd              ; data for mpi
! int_para_fc       ; integer parameters sat on each face
! nond              ; number of vertexes
! nofc              ; number of faces
! nond_on_face      ; number of vertexes which compose a face
! number_element_dof; degree of freedom on a face
! ndble_para_fc     ; number of real*8 parameters
! nint_para_fc      ; number of integer parameters
 integer function bembb2hacapk(st_bemv, st_ctl, np, face2node,dble_para_fc,rhs,sol,ztol, &
                                int_para_fc,nond,nofc,nond_on_face,number_element_dof,ndble_para_fc,nint_para_fc)
 include 'mpif.h'
 type(coordinate) np(nond)
 real*8 :: rhs(nofc*number_element_dof),sol(nofc*number_element_dof),dble_para_fc(ndble_para_fc,nofc),ztol
 integer*4 :: face2node(nond_on_face,nofc),int_para_fc(nint_para_fc,nofc)
 integer*4, dimension(:), allocatable :: lwww
 integer*4,pointer :: lpmd(:)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 real*8,dimension(:,:),allocatable :: zgmin,zgmax,zgmid
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)
 
  lpmd => st_ctl%lpmd(:)
  mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
  if(st_ctl%param(1)>1 .and. mpinr==0)  print*,'func bembb2hacapk start'
  allocate(st_bemv%int_para_fc(nint_para_fc,nofc), st_bemv%face2node(nond_on_face,nofc))
  allocate(st_bemv%dble_para_fc(ndble_para_fc,nofc))
  allocate(st_bemv%np(nond))
  st_bemv%nond=nond; st_bemv%nofc=nofc; 
  st_bemv%np=np; 
  if(nint_para_fc  /= 0) then
     st_bemv%int_para_fc=int_para_fc
  endif
     st_bemv%nint_para_fc=nint_para_fc
  if(ndble_para_fc /= 0) then
     st_bemv%dble_para_fc=dble_para_fc
  endif
  st_bemv%ndble_para_fc=ndble_para_fc; st_bemv%nond_on_face=nond_on_face; 
  st_bemv%number_element_dof=number_element_dof; st_bemv%face2node=face2node

 ndim=3
 allocate(st_bemv%zx(nond),st_bemv%zy(nond),st_bemv%zz(nond),stat=ierr)
 if(ierr/=0)then
   goto 9999
 endif 
 do il=1,nond
   st_bemv%zx(il)=np(il)%x; st_bemv%zy(il)=np(il)%y; st_bemv%zz(il)=np(il)%z; 
 enddo
 
 allocate(zgmid(nofc,ndim))
  do il=1,nofc
    n1 = face2node(1,il)+1; n2 = face2node(2,il)+1; n3 = face2node(3,il)+1;
    zgmid(il,1) = (np(n1)%x+np(n2)%x+np(n3)%x)/3.0d0; ! center of balance
    zgmid(il,2) = (np(n1)%y+np(n2)%y+np(n3)%y)/3.0d0
    zgmid(il,3) = (np(n1)%z+np(n2)%z+np(n3)%z)/3.0d0
!***************************************************
    rhs(il)=zgmid(il,3)
!***************************************************
  enddo
  bembb2hacapk= hacapk_gensolv(st_leafmtxp,st_bemv,st_ctl,zgmid,rhs,sol,ztol)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp)
  if(st_ctl%param(1)>1 .and. mpinr==0)  print*,'func bembb2hacapk end'
9999 continue
 endfunction

endmodule m_ppohBEM_bembb2hacapk
