!=====================================================================*
!                                                                     *
!   Software Name : HACApK                                            *
!         Version : 0.2.0                                             *
!                                                                     *
!   License                                                           *
!     This file is part of HACApK.                                    *
!     HACApK is a free software, you can use it under the terms       *
!     of The MIT License (MIT). See LICENSE file and User's guide     *
!     for more details.                                               *
!                                                                     *
!   ppOpen-HPC project:                                               *
!     Open Source Infrastructure for Development and Execution of     *
!     Large-Scale Scientific Applications on Post-Peta-Scale          *
!     Supercomputers with Automatic Tuning (AT).                      *
!                                                                     *
!   Sponsorship:                                                      *
!     Japan Science and Technology Agency (JST), Basic Research       *
!     Programs: CREST, Development of System Software Technologies    *
!     for post-Peta Scale High Performance Computing.                 *
!                                                                     *
!   Copyright (c) 2014 <Akihiro Ida and Takeshi Iwashita>             *
!                                                                     *
!=====================================================================*
module m_HACApK_calc_entry_ij

use m_ppohBEM_matrix_element_ij
!*** type :: st_HACApK_calc_entry
  type :: st_HACApK_calc_entry
  real*8,pointer :: ao(:)
  integer :: nd,lp61
  integer :: nond,nofc,number_element_dof,nond_on_face,nint_para_fc,ndble_para_fc,lperp
  integer,pointer :: int_para_fc(:,:), face2node(:,:)
  real*8,pointer  :: dble_para_fc(:,:)
  type(coordinate),pointer :: np(:)

  real*8,pointer :: zx(:),zy(:),zz(:)
  end type st_HACApK_calc_entry

public :: HACApK_entry_ij

contains
!***HACApK_entry_ij
  real*8 function HACApK_entry_ij(i, j, zbemv)
  type(st_HACApK_calc_entry) :: zbemv
  
  HACApK_entry_ij=ppohBEM_matrix_element_ij(i, j, &
                                            zbemv%nond, zbemv%nofc, zbemv%nond_on_face, zbemv%np, &
                                            zbemv%int_para_fc, zbemv%nint_para_fc, zbemv%dble_para_fc, &
                                            zbemv%ndble_para_fc, zbemv%face2node,zbemv%lperp)
  
  if(zbemv%lp61==3) HACApK_entry_ij=HACApK_entry_ij*zbemv%ao(i)*zbemv%ao(j)

end function HACApK_entry_ij

endmodule m_HACApK_calc_entry_ij
