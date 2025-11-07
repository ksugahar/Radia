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
