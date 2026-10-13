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
 *   The (i,j) element of the interaction matrix A = -N - diag(1/chi) (ELF-compatible)
 *
 * Note: Radia supplies problem-specific interaction entries through this
 * callback.  HACApK only sees A(i,j); the physical unknown and material
 * scaling stay on the Radia side.
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

/* Opt-in native GPU build boundary. CPU builds keep scalar evaluation.
 * Install at link time, not in mutable process-global callback state. The
 * implementation must complete synchronously, capture errors at the C ABI,
 * and never silently fall back to a different physical entry strategy. */
#ifdef RAD_HACAPK_BATCH_ENTRIES
extern void cHACApK_entry_batch(int count, const int *rows, const int *columns,
                              int i_bemv, double *values);
#endif

static inline void cHACApK_active_entry_vector(
    int ndp, int ip, int nstrtl, int nstrtt, const int *lod,
    int i_bemv, const int *mask, int direction, double *vec)
{
#ifdef RAD_HACAPK_BATCH_ENTRIES
    int rows[256], columns[256], positions[256], count=0;
    double values[256];
#endif
    int ii;
    for(ii=0;ii<ndp;++ii){
        if(mask[ii]!=0)continue;
        const int row=lod[(direction==0?ip:ii)+nstrtl];
        const int column=lod[(direction==0?ii:ip)+nstrtt];
#ifdef RAD_HACAPK_BATCH_ENTRIES
        rows[count]=row;columns[count]=column;positions[count]=ii;++count;
        if(count==256){
            int j;cHACApK_entry_batch(count,rows,columns,i_bemv,values);
            for(j=0;j<count;++j)vec[positions[j]]=values[j];count=0;
        }
#else
        vec[ii]=cHACApK_entry_ij(row,column,i_bemv);
#endif
    }
#ifdef RAD_HACAPK_BATCH_ENTRIES
    if(count){
        int j;cHACApK_entry_batch(count,rows,columns,i_bemv,values);
        for(j=0;j<count;++j)vec[positions[j]]=values[j];
    }
#endif
}

#ifdef __cplusplus
}
#endif

#endif /* CHACAPK_CALC_ENTRY_IJ_H_INCLUDED */
