/* Optional native-only dense-leaf ownership boundary.
 * Undefined by ordinary CPU and current mexcuda builds.
 * A selected consumer must copy borrowed index metadata before returning.
 * Returning nonzero takes responsibility for a still-unmaterialized leaf;
 * no matrix apply/packing/publication may occur until materialization and
 * deferred error checks succeed. C++ exceptions must never cross this ABI.
 * Zero explicitly selects the original host fill, not an error fallback.
 */
#ifndef CHACAPK_DENSE_LEAF_SINK_H
#define CHACAPK_DENSE_LEAF_SINK_H
#ifdef RAD_HACAPK_DENSE_LEAF_SINK
#ifdef __cplusplus
extern "C" {
#endif
int cHACApK_dense_leaf_sink(const st_cHACApK_leafmtx_t *leaf,
                          const int *row_order, const int *column_order,
                          int problem);
#ifdef __cplusplus
}
#endif
#endif
#endif
