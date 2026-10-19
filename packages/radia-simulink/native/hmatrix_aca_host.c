/* Existing HACApK ACA arithmetic, compiled without TaskManager/cluster/solvers.
 * CPU host compression control is intentional; entry integration runs CUDA.
 * The same included bodies remain in the production CPU HACApK base unit.
 */
#include "cHACApK_base.h"
#include "cHACApK_calc_entry_ij.h"
#include "cHACApK_lib.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cHACApK_aca.inc"
