/*
 * rad_mmm_hacapk_callback.c
 *
 * HACApK entry callback implementation for mmm_core module
 * Provides cHACApK_entry_ij function required by cHACApK_base.c
 *
 * This is a simplified version that uses a global callback function pointer
 * set from C++ code (rad_mmm_hacapk.cpp).
 */

#include <stdlib.h>

/*
 * Global callback function pointer
 * Set by C++ code before H-matrix construction
 */
typedef double (*HACApK_mmm_entry_func)(int i, int j, int bemv);
static HACApK_mmm_entry_func g_mmm_entry_func = NULL;

/**
 * Set the entry function callback from C++
 */
void HACApK_mmm_set_entry_func(HACApK_mmm_entry_func func) {
    g_mmm_entry_func = func;
}

/**
 * Clear the entry function callback
 */
void HACApK_mmm_clear_entry_func(void) {
    g_mmm_entry_func = NULL;
}

/**
 * cHACApK_entry_ij - Required by cHACApK_base.c for matrix element computation
 *
 * @param i Row index (1-based, permuted)
 * @param j Column index (1-based, permuted)
 * @param bemv BEM variable (not used in MMM)
 * @return Matrix element A(i,j)
 */
double cHACApK_entry_ij(int i, int j, int bemv) {
    if (g_mmm_entry_func != NULL) {
        return g_mmm_entry_func(i, j, bemv);
    }
    return 0.0;
}
