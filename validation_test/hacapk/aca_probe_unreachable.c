/* Link-only guards for unused tree/BLAS entry points in cHACApK_base.c.
 * The ACA probe does not call these functions; any such call must fail loud.
 * Do not link this file into a production binary.
 */
#include <stdlib.h>

void cHACApK_get_cluster_strategy(void) { abort(); }
void cHACApK_pca_split(void) { abort(); }
void hacapk_parallel_for(void) { abort(); }
void dgeqp3_(void) { abort(); }
void dorgqr_(void) { abort(); }
void dgesvd_(void) { abort(); }
