/*
!=====================================================================*
!                                                                     *
!   Software Name : HACApK                                            *
!         Version : 1.3.0                                             *
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
!   Copyright (c) 2015 <Akihiro Ida and Takeshi Iwashita>             *
!                                                                     *
!=====================================================================*
*/
#include "cHACApK_lib.h"

//***cHACApK_med3
int cHACApK_med3(
  int nl,
  int nr,
  int nlr2)
{
  if(nl < nr) {
    if (nr < nlr2) { return nr; } else if (nlr2 < nl) { return nl; } else { return nlr2; }
  } else {
    if (nlr2 < nr) { return nr;  } else if (nl < nlr2) { return nl; } else { return nlr2; }
  }
}
