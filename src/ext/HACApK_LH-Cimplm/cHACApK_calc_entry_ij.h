/*
 * HACApK Matrix Entry Callback Interface for Radia
 *
 * This header declares the callback function that HACApK uses to compute
 * matrix elements. The actual implementation is provided by Radia in
 * rad_hacapk_interface.cpp.
 *
 * The callback function signature:
 *   double cHACApK_entry_ij(int i, int j, int i_bemv)
 *
 * Parameters:
 *   i       - Global row index (1-based in original HACApK)
 *   j       - Global column index (1-based in original HACApK)
 *   i_bemv  - Problem type identifier (unused in Radia)
 *
 * Returns:
 *   The (i,j) element of the interaction matrix A = -N + diag(1/chi)
 *
 * Note: In Radia's MMM (Magnetic Moment Method), the system matrix is:
 *   A * M = H_ext + M_r/chi
 *   where A_ij = -N_ij + delta_ij/chi_i
 *         N_ij is the demagnetization tensor element
 *
 * Thread Safety:
 *   This function must be thread-safe as it may be called from multiple
 *   OpenMP threads during ACA+ compression.
 *
 * Copyright (c) 2025 Radia Project
 * License: MIT
 */

#ifndef CHACAPK_CALC_ENTRY_IJ_H_INCLUDED
#define CHACAPK_CALC_ENTRY_IJ_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Callback function to compute matrix element A(i,j)
 *
 * This function is called by HACApK during:
 *   1. ACA+ low-rank approximation
 *   2. Dense block filling
 *   3. Diagonal updates for nonlinear iteration
 *
 * The function must be implemented by the application (Radia).
 */
extern double cHACApK_entry_ij(int i, int j, int i_bemv);

#ifdef __cplusplus
}
#endif

#endif /* CHACAPK_CALC_ENTRY_IJ_H_INCLUDED */
