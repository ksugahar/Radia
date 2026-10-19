/* Optional pre-fill admission boundary. Undefined in ordinary CPU builds.
 * Leaves are borrowed, one-based, numerically empty and must not be retained
 * or mutated. Return nonzero to admit; zero rejects before any entry fill.
 * Exceptions must never cross this C ABI. The caller still owns ctl cleanup.
 */
#ifndef CHACAPK_FRAME_ADMISSION_H
#define CHACAPK_FRAME_ADMISSION_H
#ifdef RAD_HACAPK_FRAME_ADMISSION
#ifdef __cplusplus
extern "C" {
#endif
int cHACApK_frame_admission(st_cHACApK_leafmtx const *leaves,
                          int count,int dimension,int rank_cap);
#ifdef __cplusplus
}
#endif
#endif
#endif
