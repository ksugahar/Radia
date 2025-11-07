!=====================================================================*
!                                                                     *
!   Software Name : ppohBEM                                           *
!         Version : 0.1                                               *
!                                                                     *
!   License                                                           *
!     This file is part of ppohBEM.                                   *
!     ppohBEM is a free software, you can use it under the terms   *
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
!   Copyright (c) 2012 <Takeshi Iwashita, Takeshi Mifune, Yuki Noseda,*
!                    Yasuhito Takahashi, Masatoshi Kawai, Akihiro Ida>*
!                                                                     *
!=====================================================================*
!
module m_ppohBEM_matrix_element_ij
  use m_ppohBEM_user_func
  
!*** type :: coordinate
  type :: coordinate
    real(8) :: x ,y ,z
  end type coordinate
  
contains
  real(8) function ppohBEM_matrix_element_ij(i, j, nond, nofc, nond_on_fc, np,int_para_fc, nint_para_fc, dble_para_fc, ndble_para_fc, face2node, lperp)
 
  integer ,intent(in) :: i, j, nond, nofc, nond_on_fc, nint_para_fc, ndble_para_fc,lperp
  type(coordinate), intent(in) :: np(*)
  integer, intent(in) :: face2node(3, *), int_para_fc(nint_para_fc,*)
  real(8), intent(in) :: dble_para_fc(ndble_para_fc,*)
  
  integer :: n(3)
  real(8) :: xf(3), yf(3), zf(3)
  real(8) :: xp, yp, zp
  real(8) :: u(3), v(3), wi(3),wj(3) 
  real(8) :: zwidwj

  n(1:3) = face2node(1:3, i) + 1
  xf(1:3) = np( n(1:3) )%x
  yf(1:3) = np( n(1:3) )%y
  zf(1:3) = np( n(1:3) )%z

  xp = sum( xf(1:3) ) / 3d0
  yp = sum( yf(1:3) ) / 3d0
  zp = sum( zf(1:3) ) / 3d0

  if(lperp==1)then
      u(1) = xf(2) - xf(1);  v(1) = xf(3) - xf(2)
      u(2) = yf(2) - yf(1);  v(2) = yf(3) - yf(2)
      u(3) = zf(2) - zf(1);  v(3) = zf(3) - zf(2)
      call cross_product(u, v, wi)
      wi(:) = wi(:) / sqrt( dot_product(wi, wi) )
  endif

  n(1:3) = face2node(1:3, j) + 1
  xf(1:3) = np( n(1:3) )%x
  yf(1:3) = np( n(1:3) )%y
  zf(1:3) = np( n(1:3) )%z
  
  if(lperp==1)then
      u(1) = xf(2) - xf(1);  v(1) = xf(3) - xf(2)
      u(2) = yf(2) - yf(1);  v(2) = yf(3) - yf(2)
      u(3) = zf(2) - zf(1);  v(3) = zf(3) - zf(2)
      call cross_product(u, v, wj)
      wj(:) = wj(:) / sqrt( dot_product(wj, wj) )
      zwidwj=wi(1)*wj(1)+wi(2)*wj(2)+wi(3)*wj(3)
!!!  write(53,*) i,j,zwidwj
  endif

  
  if(lperp==1)then
    ppohBEM_matrix_element_ij = face_integral(xf, yf, zf, xp, yp, zp)*abs(zwidwj)
  else
    ppohBEM_matrix_element_ij = face_integral(xf, yf, zf, xp, yp, zp)
  endif
  
end function ppohBEM_matrix_element_ij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!    ppohBEM_right_hand_side_vector_element_i    !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(8) function ppohBEM_right_hand_side_vector_element_i(i, nond, nofc, nond_on_fc, np, int_para_fc, nint_para_fc, &
  dble_para_fc, ndble_para_fc, face2node)

  integer ,intent(in) :: i, nond, nofc, nond_on_fc, nint_para_fc, ndble_para_fc  !!!! call by value
  type(coordinate), intent(in) :: np(*)
  integer, intent(in) :: face2node(3, *), int_para_fc(nint_para_fc,*)
  real(8), intent(in) :: dble_para_fc( ndble_para_fc, * )
  
  ppohBEM_right_hand_side_vector_element_i = dble_para_fc(1,i)
  
end function ppohBEM_right_hand_side_vector_element_i

end module m_ppohBEM_matrix_element_ij
